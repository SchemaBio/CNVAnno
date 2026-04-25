"""
CNVAnno - ClinGen规则CNV致病性评分注释器

基于ClinGen CNV Working Group官方评分规则进行CNV致病性评估。
参考: https://cnvcalc.clinicalgenome.org/

评分系统采用累积评分制:
- 总分 ≤ -0.99: Benign (良性)
- 总分 -0.98 ~ 0.99: VUS (意义未明)
- 总分 ≥ 1.00: Pathogenic (致病)
- 总分 ≥ 1.99: Likely Pathogenic (可能致病，有些版本)

Author: CNVAnno
Date: 2025-04-25
"""

import argparse
import gzip
import json
import logging
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd
from intervaltree import IntervalTree

# ============================================================================
# ClinGen CNV评分规则常量
# ============================================================================

# 最终分类阈值 (ClinGen官方)
CLASSIFICATION_THRESHOLDS = {
    "Benign": -0.99,           # 总分 ≤ -0.99
    "Likely Benign": -0.98,    # 有些版本使用 -0.98 作为Likely Benign阈值
    "VUS_lower": -0.98,        # VUS区间下界
    "VUS_upper": 0.99,         # VUS区间上界 (-0.98 ~ 0.99)
    "Likely Pathogenic": 1.00, # 有些版本使用 1.00 作为Likely Pathogenic阈值
    "Pathogenic": 1.00,        # 总分 ≥ 1.00
}

# Section 1: 基因组内容初始评估
SECTION1_SCORES = {
    "contains_protein_coding": 0.0,      # 1A: 含蛋白编码基因
    "no_protein_coding": -0.60,          # 1B: 不含蛋白编码基因
}

# Section 2: 已知基因/区域重叠评分
# Loss (Deletion) 评分
SECTION2_LOSS_SCORES = {
    "2A_full_overlap_HI": 1.0,           # 完全重叠已确立HI基因/区域
    "2B_partial_overlap_HI_region": 0.0, # 部分重叠HI区域(未涉及致病基因)
    "2C_1_partial_5prime_coding": 0.90,  # 部分重叠5'端，涉及编码序列 (最大)
    "2C_2_partial_5prime_utr": 0.45,     # 部分重叠5'端，仅5' UTR
    "2D_2_partial_3prime_last_exon_pathogenic": 0.90, # 部分3'端，最后外显子(有致病报道)
    "2D_3_partial_3prime_last_exon_no_pathogenic": 0.45, # 部分3'端，最后外显子(无致病报道)
    "2D_4_partial_3prime_multiple_exons": 0.90, # 部分3'端，含多个外显子
    "2E_both_breakpoints_in_gene": 0.90, # 两断点在同一基因内(PVS1)
    "2F_in_benign_region": -1.0,         # 完全位于良性CNV区域
    "2H_multiple_HI_predictors": 0.15,   # 多个HI预测器提示
}

# Gain (Duplication) 评分
SECTION2_GAIN_SCORES = {
    "2A_full_overlap_TS": 1.0,           # 完全重叠已确立TS基因/区域
    "2B_partial_overlap_TS": 0.0,        # 部分重叠TS区域
    "2C_same_genes_as_benign": -1.0,     # 基因内容与良性CNV相同
    "2D_smaller_no_gene_disruption": -1.0, # 较小，不中断基因
    "2E_smaller_possible_disruption": 0.0, # 较小，可能中断基因
    "2F_larger_no_extra_genes": -0.90,   # 较大，无额外基因
    "2G_larger_with_extra_genes": 0.0,   # 较大，含额外基因
    "2I_both_breakpoints_in_HI_gene_PVS1": 0.90, # 两断点在HI基因内(PVS1)
    "2K_one_breakpoint_HI_gene_phenotype": 0.45, # 一断点在HI基因+表型一致
}

# Section 3: 基因数量评分
SECTION3_SCORES = {
    "loss": {
        0: 0.0,      # 0基因
        1: 0.0,      # 1基因
        2: 0.15,     # 2-4基因
        5: 0.30,     # 5-9基因
        10: 0.45,    # 10-24基因
        25: 0.60,    # 25-49基因
        50: 0.90,    # ≥50基因
    },
    "gain": {
        0: 0.0,      # 0-34基因
        35: 0.45,    # 35-49基因
        50: 0.90,    # ≥50基因
    },
}

# Section 4: 文献和数据库证据
# 注意: 以下证据需要人工决策(表型、新发、分离等)，标记为"MANUAL"
SECTION4_SCORES = {
    # De novo证据 - 需要人工决策
    "4A_de_novo_highly_specific": "MANUAL: Confirm de novo + highly specific phenotype (0.45/assumed 0.30)",
    "4B_de_novo_consistent_specific": "MANUAL: Confirm de novo + consistent phenotype (0.30/assumed 0.15)",
    "4C_de_novo_consistent_not_specific": "MANUAL: Confirm de novo + non-specific phenotype (0.15/assumed 0.10)",
    "4D_de_novo_inconsistent": "MANUAL: De novo with inconsistent phenotype (0 to -0.30)",

    # 分离分析证据 - 需要人工决策
    "4F_segregation": "MANUAL: Co-segregates with phenotype (0.15-0.45 depending on occurrences)",
    "4G_segregation_2_occurrences": 0.15,
    "4H_segregation_3_occurrences": 0.30,

    # 非分离分析 - 需要人工决策
    "4I_affected_not_carrying": "MANUAL: Affected relative does not carry (-0.45)",
    "4J_unaffected_carrying_specific": "MANUAL: Unaffected parent carries (specific phenotype) (-0.30)",
    "4K_unaffected_carrying_not_specific": "MANUAL: Unaffected parent carries (non-specific phenotype) (-0.15)",

    # 病例对照研究 - 需要人工决策
    "4L_case_control_significant": "MANUAL: Case-control significant (specific phenotype) (0.45)",
    "4M_case_control_not_specific": "MANUAL: Case-control significant (non-specific) (0.30)",
    "4N_case_control_no_difference": "MANUAL: Case-control no significant difference (-0.90)",

    # 频率证据 - 可自动评分
    "4O_overlap_with_common_population": -1.00,  # 与常见人群变异重叠 (>1%)
    "4O_moderate_frequency": -0.45,              # 中等频率 (0.1%-1%)
}

# Section 5: 遗传模式评分 - 需要人工决策
SECTION5_SCORES = {
    "5A_de_novo": "MANUAL: Use Section 4 de novo scores",
    "5B_inherited_unaffected_specific": "MANUAL: Inherited from unaffected parent (specific phenotype) (0 to -0.45)",
    "5C_inherited_unaffected_not_specific": "MANUAL: Inherited from unaffected parent (non-specific) (0 to -0.30)",
    "5D_segregates": "MANUAL: Segregates with family phenotype (use Section 4F-H)",
    "5F_no_info": 0.0,  # 无遗传信息
    "5G_no_info_not_specific": "MANUAL: No info available, non-specific phenotype (0-0.15)",
    "5H_no_info_highly_specific": "MANUAL: No info available, highly specific phenotype (0-0.30)",
}

# 频率阈值定义
FREQUENCY_THRESHOLDS = {
    "benign": 0.01,           # >1% 触发4O证据(score=-1.00)
    "likely_benign": 0.001,   # >0.1% 部分良性证据(score=-0.45)
}

# ClinGen HI/TS评分含义映射
HI_SCORE_MEANINGS = {
    3: "Sufficient evidence for haploinsufficiency (established HI gene)",
    2: "Likely haploinsufficient",
    1: "Little evidence for haploinsufficiency",
    0: "No evidence available",
    30: "No evidence available (gene not curated)",
    40: "Dosage sensitivity unlikely (benign evidence)",
}

TS_SCORE_MEANINGS = {
    3: "Sufficient evidence for triplosensitivity (established TS gene)",
    2: "Likely triplosensitive",
    1: "Little evidence for triplosensitivity",
    0: "No evidence available",
    40: "Dosage sensitivity unlikely (benign evidence)",
}

# 外显子信息结构
@dataclass
class ExonInfo:
    """外显子信息"""
    gene: str
    transcript: str
    exon_number: int
    strand: str  # '+' or '-'
    chrom: str
    start: int
    end: int
    cytoband: str

    def is_first_exon(self) -> bool:
        """是否是第一外显子(5'端)"""
        if self.strand == '+':
            return self.exon_number == 1
        else:  # '-' strand: exon编号大的是5'端
            # 需要知道总外显子数才能判断
            return False  # 需要额外信息

    def get_5prime_position(self) -> str:
        """获取5'端位置描述"""
        if self.strand == '+':
            return f"exon {self.exon_number} (5' end for + strand)"
        else:
            return f"exon {self.exon_number} (3' end for - strand)"


# ============================================================================
# 数据类定义
# ============================================================================

def normalize_chrom(chrom: str) -> str:
    """标准化染色体名称为 chrN 格式"""
    chrom = str(chrom).strip()
    chrom_lower = chrom.lower()
    if chrom_lower.startswith('chr'):
        chrom_num = chrom[3:]
    else:
        chrom_num = chrom
    chrom_num_lower = chrom_num.lower()
    if chrom_num_lower in ['mt', 'm']:
        return 'chrM'
    try:
        int(chrom_num)
        return f'chr{chrom_num}'
    except ValueError:
        return f'chr{chrom_num.upper()}'


@dataclass
class CNVRecord:
    """单个CNV记录"""
    chrom: str
    start: int
    end: int
    status: str  # Deletion/Amplification/Normal
    size: int = 0

    def __post_init__(self):
        self.size = self.end - self.start
        self.chrom = normalize_chrom(self.chrom)


@dataclass
class GeneDosageInfo:
    """基因剂量敏感性信息"""
    gene: str
    hi_score: int
    hi_desc: str
    tr_score: int
    tr_desc: str
    chrom: str
    start: int
    end: int


@dataclass
class RegionCurationInfo:
    """区域级别的剂量敏感性信息"""
    isca_id: str
    region_name: str
    cytoband: str
    chrom: str
    start: int
    end: int
    hi_score: int
    hi_desc: str
    tr_score: int
    tr_desc: str


@dataclass
class FrequencyInfo:
    """人群频率信息"""
    source: str
    deletion_freq: float
    duplication_freq: float
    sample_size: int


@dataclass
class GenCCInfo:
    """基因-疾病关联信息"""
    gene: str
    disease: str
    classification: str
    moi: str


@dataclass
class ClinGenEvidence:
    """ClinGen证据评分详情"""
    # 各Section分数
    section1_score: float = 0.0
    section2_score: float = 0.0
    section3_score: float = 0.0
    section4_score: float = 0.0
    section5_score: float = 0.0
    total_score: float = 0.0

    # 客观证据使用标记 (可自动评分)
    evidence_1A: int = 0  # 含蛋白编码基因
    evidence_1B: int = 0  # 不含蛋白编码基因
    evidence_2A: int = 0  # 完全重叠致病基因/区域 (HI/TR=3)
    evidence_2B: int = 0  # 部分重叠致病基因
    evidence_2C: int = 0  # 部分5'端重叠
    evidence_2D: int = 0  # 部分3'端重叠
    evidence_2E: int = 0  # 两断点在同一基因内
    evidence_2F: int = 0  # 位于良性区域(HI/TR=40)
    evidence_2H: int = 0  # Likely HI/TS基因(score=2)
    evidence_2K: int = 0  # Gain: 一断点在HI基因内(需表型确认)
    evidence_3_gene_count: int = 0  # 基因数量评分
    evidence_4O: int = 0  # 高人群频率(良性)
    evidence_4A: int = 0  # de novo证据(需要人工决策)
    evidence_4L: int = 0  # 病例对照研究(需要人工决策)
    evidence_5: int = 0  # 遗传模式(需要人工决策)

    # 主观证据标记 (需要人工决策)
    manual_decision_needed: int = 0  # 是否需要人工决策
    manual_4A_de_novo: str = ""  # 证据4A: 新发变异(需确认+表型特异性)
    manual_4B_de_novo: str = ""  # 证据4B: 新发+表型一致
    manual_4C_de_novo: str = ""  # 证据4C: 新发+表型非特异
    manual_4D_de_novo: str = ""  # 证据4D: 新发+表型不一致(负分)
    manual_4F_segregation: str = ""  # 证据4F-H: 分离分析
    manual_4I_non_segregation: str = ""  # 证据4I-K: 非分离(负分)
    manual_4L_case_control: str = ""  # 证据4L-M: 病例对照研究
    manual_5B_inheritance: str = ""  # 证据5B-C: 遗传自未受累父母(负分)
    manual_5D_segregation: str = ""  # 证据5D: 家系分离
    manual_5H_specific_phenotype: str = ""  # 证据5H: 无遗传信息+特异表型

    # 详细证据描述
    evidence_details: List[str] = field(default_factory=list)

    # 外显子重叠详情 (用于Section 2C/2D)
    exon_overlap_details: List[str] = field(default_factory=list)


@dataclass
class CNVAnnotation:
    """CNV完整注释结果"""
    cnv: CNVRecord
    overlapping_genes: List[str] = field(default_factory=list)
    gene_count: int = 0
    dosage_sensitive_genes: List[GeneDosageInfo] = field(default_factory=list)
    overlapping_regions: List[RegionCurationInfo] = field(default_factory=list)
    frequency_info: List[FrequencyInfo] = field(default_factory=list)
    gencc_info: List[GenCCInfo] = field(default_factory=list)
    hi_max_score: int = 0
    tr_max_score: int = 0
    max_freq: float = 0.0
    gencc_supportive_genes: List[str] = field(default_factory=list)

    # 外显子重叠信息 (用于Section 2C/2D评分)
    overlapping_exons: List[ExonInfo] = field(default_factory=list)
    affected_first_exons: List[str] = field(default_factory=list)  # 受影响的第一外显子基因
    affected_last_exons: List[str] = field(default_factory=list)   # 受影响的最后外显子基因

    # ClinGen评分
    clingen_evidence: ClinGenEvidence = field(default_factory=ClinGenEvidence)
    overall_classification: str = "Uncertain Significance"
    classification_reason: str = ""

    # ISCN标准命名
    iscn: str = ""


# ============================================================================
# ISCN格式化器
# ============================================================================

class ISCNFormatter:
    """生成ISCN标准命名格式"""

    def __init__(self, genome_build: str = "GRCh38"):
        self.genome_build = genome_build

    def format_iscn(self, chrom: str, start: int, end: int, status: str,
                    start_cytoband: str = "", end_cytoband: str = "") -> str:
        """
        生成ISCN标准命名

        Args:
            chrom: 染色体号 (如 "1", "chr1", "X")
            start: 起始坐标
            end: 终止坐标
            status: CNV状态 (Deletion/Amplification等)
            start_cytoband: 起始位置的cytoband
            end_cytoband: 终止位置的cytoband

        Returns:
            ISCN格式字符串
        """
        # 标准化染色体号
        chrom_clean = chrom.replace("chr", "")

        # 确定CNV类型
        if status.lower() in ["deletion", "del", "loss"]:
            cnv_type = "del"
        elif status.lower() in ["amplification", "dup", "duplication", "gain"]:
            cnv_type = "dup"
        else:
            cnv_type = "cnv"  # 其他类型用通用cnv表示

        # 构建ISCN
        # 格式: seq[GRCh38] 染色体号(起始_终止)del 或 带号写法
        build_name = self.genome_build

        if start_cytoband and end_cytoband:
            # 有cytoband信息，使用带号写法
            # 简化带号：如果起止在同一带，只写一个带号
            if start_cytoband == end_cytoband:
                band_part = f"{chrom_clean}{start_cytoband}"
            else:
                band_part = f"{chrom_clean}{start_cytoband}{end_cytoband}"
            iscn = f"seq[{build_name}] {band_part}({start}_{end}){cnv_type}"
        else:
            # 无cytoband信息，使用简写格式
            iscn = f"seq[{build_name}] {chrom_clean}({start}_{end}){cnv_type}"

        return iscn


# ============================================================================
# 数据库加载器
# ============================================================================

class DatabaseLoader:
    """加载和索引各个数据库"""

    def __init__(self, data_dir: str, genome_build: str = "GRCh38"):
        self.data_dir = Path(data_dir)
        self.genome_build = genome_build
        self.hi_genes: Dict[str, GeneDosageInfo] = {}
        self.tr_genes: Dict[str, GeneDosageInfo] = {}
        self.curated_regions: List[RegionCurationInfo] = []
        self.gencc_data: Dict[str, List[GenCCInfo]] = defaultdict(list)
        self.benign_regions: List[RegionCurationInfo] = []

        # 外显子数据
        self.exon_data: Dict[str, List[ExonInfo]] = defaultdict(list)  # 按基因存储
        self.gene_exon_count: Dict[str, int] = {}  # 每个基因的外显子总数

        # Cytoband数据 (用于ISCN格式化)
        self.cytoband_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)

        self.gene_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        self.region_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        self.benign_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        self.exon_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)  # 外显子区间树
        self.gnomad_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        self.population_interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)

        self.logger = logging.getLogger(__name__)

    def load_all(self):
        """加载所有数据库"""
        self.logger.info("Loading ClinGen databases...")
        self._load_haploinsufficiency_genes()
        self._load_triplosensitivity_genes()
        self._build_gene_interval_trees()
        self._load_curated_regions()
        self._build_region_interval_trees()
        self._load_gencc()
        self._load_cytoband()  # 加载cytoband数据（用于ISCN格式化）
        self._load_gencode_exons()  # 加载Gencode外显子数据
        self._build_exon_interval_trees()
        self._load_gnomad_cnv()
        self._build_gnomad_interval_trees()
        self._load_population_cnv()
        self._build_population_interval_trees()
        self.logger.info("All databases loaded successfully")

    def _load_haploinsufficiency_genes(self):
        """加载单倍体不足基因数据库"""
        hi_file = self.data_dir / f"ClinGen_haploinsufficiency_gene_{self.genome_build}.bed"
        if not hi_file.exists():
            self.logger.warning(f"HI gene file not found: {hi_file}")
            return

        with open(hi_file) as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = normalize_chrom(parts[0])
                    start = int(parts[1])
                    end = int(parts[2])
                    gene = parts[3]
                    score = int(parts[4])
                    hi_desc = HI_SCORE_MEANINGS.get(score, "Unknown")
                    self.hi_genes[gene] = GeneDosageInfo(
                        gene=gene,
                        hi_score=score,
                        hi_desc=hi_desc,
                        tr_score=0,
                        tr_desc="",
                        chrom=chrom,
                        start=start,
                        end=end
                    )
        self.logger.info(f"Loaded {len(self.hi_genes)} haploinsufficiency genes")

    def _load_triplosensitivity_genes(self):
        """加载三倍体敏感性基因数据库"""
        tr_file = self.data_dir / f"ClinGen_triplosensitivity_gene_{self.genome_build}.bed"
        if not tr_file.exists():
            self.logger.warning(f"TR gene file not found: {tr_file}")
            return

        with open(tr_file) as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = normalize_chrom(parts[0])
                    start = int(parts[1])
                    end = int(parts[2])
                    gene = parts[3]
                    score = int(parts[4])
                    tr_desc = TS_SCORE_MEANINGS.get(score, "Unknown")

                    if gene in self.hi_genes:
                        self.hi_genes[gene].tr_score = score
                        self.hi_genes[gene].tr_desc = tr_desc
                    else:
                        self.tr_genes[gene] = GeneDosageInfo(
                            gene=gene,
                            hi_score=0,
                            hi_desc="",
                            tr_score=score,
                            tr_desc=tr_desc,
                            chrom=chrom,
                            start=start,
                            end=end
                        )

        for gene, info in self.tr_genes.items():
            if gene not in self.hi_genes:
                self.hi_genes[gene] = info

        self.logger.info(f"Total genes with dosage info: {len(self.hi_genes)}")

    def _build_gene_interval_trees(self):
        """构建基因区间索引树"""
        for gene, info in self.hi_genes.items():
            self.gene_interval_trees[info.chrom].addi(info.start, info.end, info)
        self.logger.info("Built gene interval trees")

    def _load_curated_regions(self):
        """加载ClinGen已注释的CNV区域"""
        region_file = self.data_dir / f"ClinGen_region_curation_list_{self.genome_build}.tsv"
        if not region_file.exists():
            self.logger.warning(f"Region file not found: {region_file}")
            return

        with open(region_file) as f:
            for line in f:
                if line.startswith('#ClinGen') or line.startswith('#Genomic'):
                    continue
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                genomic_loc = parts[3]
                if genomic_loc == 'tbd' or not genomic_loc:
                    continue

                match = re.match(r'chr(\w+):(\d+)-(\d+)', genomic_loc)
                if not match:
                    continue

                chrom = normalize_chrom(f'chr{match.group(1)}')
                start = int(match.group(2))
                end = int(match.group(3))

                hi_score = int(parts[4]) if parts[4] else 0
                hi_desc = parts[5] if len(parts) > 5 else ""
                tr_score = int(parts[12]) if len(parts) > 12 and parts[12] else 0
                tr_desc = parts[13] if len(parts) > 13 else ""

                region = RegionCurationInfo(
                    isca_id=parts[0],
                    region_name=parts[1],
                    cytoband=parts[2],
                    chrom=chrom,
                    start=start,
                    end=end,
                    hi_score=hi_score,
                    hi_desc=hi_desc,
                    tr_score=tr_score,
                    tr_desc=tr_desc
                )

                # 区分致病区域和良性区域
                if hi_score == 40 or tr_score == 40:
                    self.benign_regions.append(region)
                else:
                    self.curated_regions.append(region)

        self.logger.info(f"Loaded {len(self.curated_regions)} curated pathogenic regions")
        self.logger.info(f"Loaded {len(self.benign_regions)} benign regions")

    def _build_region_interval_trees(self):
        """构建区域区间索引树"""
        for region in self.curated_regions:
            self.region_interval_trees[region.chrom].addi(region.start, region.end, region)
        for region in self.benign_regions:
            self.benign_interval_trees[region.chrom].addi(region.start, region.end, region)
        self.logger.info("Built region interval trees")

    def _load_gencc(self):
        """加载GenCC基因-疾病关联数据"""
        gencc_file = self.data_dir / "gencc-submissions.xlsx"
        if not gencc_file.exists():
            self.logger.warning(f"GenCC file not found: {gencc_file}")
            return

        df = pd.read_excel(gencc_file)
        valid_classifications = ['Definitive', 'Strong', 'Moderate', 'Limited', 'Supportive']
        df = df[df['classification_title'].isin(valid_classifications)]

        for _, row in df.iterrows():
            gene = row['gene_symbol']
            self.gencc_data[gene].append(GenCCInfo(
                gene=gene,
                disease=row['disease_title'],
                classification=row['classification_title'],
                moi=row['moi_title']
            ))

        self.logger.info(f"Loaded GenCC data for {len(self.gencc_data)} genes")

    def _load_cytoband(self):
        """加载UCSC cytoband数据（用于ISCN格式化）"""
        cytoband_file = self.data_dir / f"cytoBand_{self.genome_build}.txt"
        if not cytoband_file.exists():
            self.logger.warning(f"Cytoband file not found: {cytoband_file}")
            return

        count = 0
        with open(cytoband_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue

                # 格式: chrom, start, end, band_name, gieStain
                chrom = normalize_chrom(parts[0])
                start = int(parts[1])
                end = int(parts[2])
                band_name = parts[3]  # 如 "p36.33"

                # 构建带号：将染色体号和带名组合，如 "p36.33" -> 用于ISCN
                # 注意：band_name已经包含p/q信息，ISCN格式需要组合
                # 例如：chr1 + p36.33 -> 在ISCN中显示为 "1p36.33"
                self.cytoband_interval_trees[chrom].addi(start, end, band_name)
                count += 1

        self.logger.info(f"Loaded {count} cytoband records")

    def _load_gencode_exons(self):
        """加载Gencode外显子坐标数据"""
        gencode_file = self.data_dir / f"Gencode.{self.genome_build}.cnvkit.target.bed"
        if not gencode_file.exists():
            self.logger.warning(f"Gencode exon file not found: {gencode_file}")
            return

        count = 0
        gene_exon_max: Dict[str, int] = {}  # 记录每个基因的最大外显子编号

        with open(gencode_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue

                # 格式: chrom, start, end, GENE|Transcript|ENST|ExonNum|Strand|Cytoband
                chrom = normalize_chrom(parts[0])
                start = int(parts[1])
                end = int(parts[2])

                # 解析基因信息
                gene_info = parts[3].split('|')
                if len(gene_info) < 5:
                    continue

                gene = gene_info[0]
                transcript = gene_info[1] if len(gene_info) > 1 else ""
                exon_num = int(gene_info[3]) if len(gene_info) > 3 else 0
                strand = gene_info[4] if len(gene_info) > 4 else '+'
                cytoband = gene_info[5] if len(gene_info) > 5 else ""

                exon_info = ExonInfo(
                    gene=gene,
                    transcript=transcript,
                    exon_number=exon_num,
                    strand=strand,
                    chrom=chrom,
                    start=start,
                    end=end,
                    cytoband=cytoband
                )

                self.exon_data[gene].append(exon_info)

                # 记录最大外显子编号
                if gene not in gene_exon_max or exon_num > gene_exon_max[gene]:
                    gene_exon_max[gene] = exon_num

                count += 1

        # 设置每个基因的外显子总数
        for gene, max_exon in gene_exon_max.items():
            self.gene_exon_count[gene] = max_exon

        self.logger.info(f"Loaded {count} exon records for {len(self.exon_data)} genes")

    def _build_exon_interval_trees(self):
        """构建外显子区间索引树"""
        for gene, exons in self.exon_data.items():
            for exon in exons:
                self.exon_interval_trees[exon.chrom].addi(exon.start, exon.end, exon)
        self.logger.info("Built exon interval trees")

    def get_cytoband(self, chrom: str, position: int) -> str:
        """
        根据坐标获取cytoband名称

        Args:
            chrom: 染色体号 (如 "chr1", "1")
            position: 基因组坐标

        Returns:
            cytoband名称 (如 "p36.33")，未找到返回空字符串
        """
        chrom = normalize_chrom(chrom)
        tree = self.cytoband_interval_trees.get(chrom)
        if not tree:
            return ""

        # 查找重叠的区间 (使用overlap方法)
        overlaps = tree.overlap(position, position + 1)
        if overlaps:
            # 返回第一个找到的cytoband
            for interval in overlaps:
                return interval.data

        return ""

    def _load_gnomad_cnv(self):
        """加载gnomAD CNV频率数据"""
        gnomad_file = self.data_dir / f"gnomad.v4.1.cnv.all.{self.genome_build.lower()}.vcf.gz"
        if not gnomad_file.exists():
            self.logger.warning(f"gnomAD CNV file not found: {gnomad_file}")
            return

        count = 0
        with gzip.open(gnomad_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue

                chrom = normalize_chrom(parts[0])
                pos = int(parts[1])
                info = parts[7]

                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, val = item.split('=')
                        info_dict[key] = val

                end = int(info_dict.get('END', pos))
                svtype = info_dict.get('SVTYPE', '')
                sf = float(info_dict.get('SF', 0))

                freq_info = {
                    'source': 'gnomAD',
                    'svtype': svtype,
                    'freq': sf,
                }
                self.gnomad_interval_trees[chrom].addi(pos, end, freq_info)
                count += 1

        self.logger.info(f"Loaded {count} gnomAD CNV records")

    def _build_gnomad_interval_trees(self):
        pass

    def _load_population_cnv(self):
        """加载DECIPHER人群CNV数据"""
        pop_file = self.data_dir / f"population_cnv_{self.genome_build.lower()}.txt.gz"
        if not pop_file.exists():
            self.logger.warning(f"Population CNV file not found: {pop_file}")
            return

        count = 0
        with gzip.open(pop_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue

                chrom = normalize_chrom(f'chr{parts[1]}')
                start = int(parts[2])
                end = int(parts[3])
                del_freq = float(parts[5]) if parts[5] else 0
                dup_freq = float(parts[8]) if parts[8] else 0
                sample_size = int(parts[14]) if len(parts) > 14 else 0

                freq_info = {
                    'source': 'DECIPHER',
                    'del_freq': del_freq,
                    'dup_freq': dup_freq,
                    'sample_size': sample_size
                }
                self.population_interval_trees[chrom].addi(start, end, freq_info)
                count += 1

        self.logger.info(f"Loaded {count} population CNV records")

    def _build_population_interval_trees(self):
        pass

    def query_genes_by_interval(self, chrom: str, start: int, end: int) -> List[GeneDosageInfo]:
        """通过区间查询重叠的基因"""
        results = []
        tree = self.gene_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)
        return results

    def query_regions_by_interval(self, chrom: str, start: int, end: int) -> List[RegionCurationInfo]:
        """查询重叠的致病区域"""
        results = []
        tree = self.region_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)
        return results

    def query_benign_regions_by_interval(self, chrom: str, start: int, end: int) -> List[RegionCurationInfo]:
        """查询重叠的良性区域"""
        results = []
        tree = self.benign_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)
        return results

    def query_frequency_by_interval(self, chrom: str, start: int, end: int, cnv_type: str) -> List[dict]:
        """查询人群频率"""
        results = []

        # gnomAD
        tree = self.gnomad_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)

        # DECIPHER
        tree = self.population_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)

        return results

    def query_exons_by_interval(self, chrom: str, start: int, end: int) -> List[ExonInfo]:
        """查询重叠的外显子"""
        results = []
        tree = self.exon_interval_trees.get(chrom)
        if tree:
            intervals = tree.overlap(start, end)
            for interval in intervals:
                results.append(interval.data)
        return results

    def get_gene_exon_count(self, gene: str) -> int:
        """获取基因的外显子总数"""
        return self.gene_exon_count.get(gene, 0)

    def is_first_exon(self, exon: ExonInfo) -> bool:
        """判断是否是第一外显子(转录起始端)"""
        total_exons = self.get_gene_exon_count(exon.gene)
        if exon.strand == '+':
            # +链: exon 1是5'端(转录起始)
            return exon.exon_number == 1
        else:
            # -链: exon编号最大的是5'端
            return exon.exon_number == total_exons

    def is_last_exon(self, exon: ExonInfo) -> bool:
        """判断是否是最后外显子(转录终止端)"""
        total_exons = self.get_gene_exon_count(exon.gene)
        if exon.strand == '+':
            # +链: 编号最大的是3'端(转录终止)
            return exon.exon_number == total_exons
        else:
            # -链: exon 1是3'端
            return exon.exon_number == 1


# ============================================================================
# ClinGen评分器
# ============================================================================

class ClinGenScorer:
    """基于ClinGen官方规则的CNV评分器"""

    def __init__(self, db_loader: DatabaseLoader):
        self.db = db_loader
        self.logger = logging.getLogger(__name__)

    def score_cnv(self, annotation: CNVAnnotation) -> ClinGenEvidence:
        """按ClinGen官方规则评分"""
        cnv = annotation.cnv
        evidence = ClinGenEvidence()

        if cnv.status == "Normal":
            evidence.evidence_details.append("Normal status - no CNV")
            annotation.overall_classification = "Benign"
            annotation.classification_reason = "No CNV detected"
            return evidence

        # Section 1: 基因组内容评估
        self._score_section1(annotation, evidence)

        # Section 2: 基因/区域重叠评估
        self._score_section2(annotation, evidence)

        # Section 3: 基因数量评估
        self._score_section3(annotation, evidence)

        # Section 4: 频率/数据库证据
        self._score_section4(annotation, evidence)

        # Section 5: 遗传模式 (数据不可用，标记为0)
        evidence.section5_score = 0.0
        evidence.evidence_5 = 0
        evidence.evidence_details.append("Section 5: No inheritance data available (score=0)")

        # 计算总分
        evidence.total_score = (
            evidence.section1_score +
            evidence.section2_score +
            evidence.section3_score +
            evidence.section4_score +
            evidence.section5_score
        )

        # 最终分类
        self._classify(annotation, evidence)

        return evidence

    def _score_section1(self, annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """Section 1: 基因组内容评估"""
        cnv = annotation.cnv

        # 查询重叠基因
        overlapping_genes = self.db.query_genes_by_interval(cnv.chrom, cnv.start, cnv.end)
        annotation.overlapping_genes = [g.gene for g in overlapping_genes]
        annotation.gene_count = len(annotation.overlapping_genes)
        annotation.dosage_sensitive_genes = overlapping_genes

        if annotation.gene_count > 0:
            evidence.section1_score = SECTION1_SCORES["contains_protein_coding"]
            evidence.evidence_1A = 1
            evidence.evidence_details.append(f"Section 1A: Contains {annotation.gene_count} protein-coding genes (score=0)")
        else:
            evidence.section1_score = SECTION1_SCORES["no_protein_coding"]
            evidence.evidence_1B = 1
            evidence.evidence_details.append(f"Section 1B: No protein-coding genes (score=-0.60)")

    def _score_section2(self, annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """Section 2: 基因/区域重叠评估"""
        cnv = annotation.cnv
        overlapping_genes = annotation.dosage_sensitive_genes  # 从Section1获取

        # 查询重叠的外显子
        overlapping_exons = self.db.query_exons_by_interval(cnv.chrom, cnv.start, cnv.end)
        annotation.overlapping_exons = overlapping_exons

        # 分析外显子重叠情况
        for exon in overlapping_exons:
            total_exons = self.db.get_gene_exon_count(exon.gene)
            if self.db.is_first_exon(exon):
                annotation.affected_first_exons.append(exon.gene)
                evidence.exon_overlap_details.append(
                    f"{exon.gene}: first exon (exon {exon.exon_number}, strand {exon.strand})"
                )
            if self.db.is_last_exon(exon):
                annotation.affected_last_exons.append(exon.gene)
                evidence.exon_overlap_details.append(
                    f"{exon.gene}: last exon (exon {exon.exon_number}, strand {exon.strand})"
                )

        # 检查良性区域重叠 (Section 2F)
        benign_regions = self.db.query_benign_regions_by_interval(cnv.chrom, cnv.start, cnv.end)
        for region in benign_regions:
            overlap_pct = self._calculate_overlap_pct(cnv, region)
            if overlap_pct >= 0.80:  # 80%以上重叠良性区域
                evidence.section2_score = SECTION2_LOSS_SCORES["2F_in_benign_region"]
                evidence.evidence_2F = 1
                evidence.evidence_details.append(f"Section 2F: Overlaps benign region {region.region_name} (score=-1.00)")
                annotation.overlapping_regions.append(region)
                return  # 良性证据优先，返回

        # 检查致病区域重叠
        pathogenic_regions = self.db.query_regions_by_interval(cnv.chrom, cnv.start, cnv.end)
        annotation.overlapping_regions = pathogenic_regions

        # 根据CNV类型评分
        if cnv.status == "Deletion":
            self._score_section2_loss(annotation, evidence, pathogenic_regions, overlapping_genes)
        elif cnv.status == "Amplification":
            self._score_section2_gain(annotation, evidence, pathogenic_regions, overlapping_genes)

    def _score_section2_loss(self, annotation: CNVAnnotation, evidence: ClinGenEvidence,
                             pathogenic_regions: List, overlapping_genes: List):
        """Loss (Deletion) Section 2评分 - 使用外显子信息"""
        cnv = annotation.cnv

        # 检查完全重叠HI score=3的基因
        for gene_info in overlapping_genes:
            if gene_info.hi_score == 3:
                # 检查是否完全重叠
                if cnv.start <= gene_info.start and cnv.end >= gene_info.end:
                    evidence.section2_score = SECTION2_LOSS_SCORES["2A_full_overlap_HI"]
                    evidence.evidence_2A = 1
                    evidence.evidence_details.append(
                        f"Section 2A: Full overlap with established HI gene {gene_info.gene} (score=1.00)"
                    )
                    annotation.hi_max_score = 3
                    return
                else:
                    # 部分重叠 - 需要分析外显子情况
                    self._analyze_partial_overlap_loss(gene_info, annotation, evidence)
                    annotation.hi_max_score = max(annotation.hi_max_score, gene_info.hi_score)

        # 检查致病区域重叠
        for region in pathogenic_regions:
            if region.hi_score == 3:
                overlap_pct = self._calculate_overlap_pct(cnv, region)
                if overlap_pct >= 0.80:
                    evidence.section2_score = SECTION2_LOSS_SCORES["2A_full_overlap_HI"]
                    evidence.evidence_2A = 1
                    evidence.evidence_details.append(
                        f"Section 2A: High overlap with pathogenic region {region.region_name} (score=1.00)"
                    )
                    annotation.hi_max_score = 3
                    return
                elif overlap_pct >= 0.50:
                    evidence.section2_score = max(evidence.section2_score, 0.45)
                    evidence.evidence_2B = 1
                    evidence.evidence_details.append(
                        f"Section 2B: Partial overlap with region {region.region_name} (score=0.45)"
                    )

        # HI score=2 的基因
        for gene_info in overlapping_genes:
            if gene_info.hi_score == 2:
                evidence.section2_score = max(evidence.section2_score, 0.15)
                evidence.evidence_2H = 1
                evidence.evidence_details.append(
                    f"Section 2H: Likely HI gene {gene_info.gene} (score=0.15)"
                )
                annotation.hi_max_score = max(annotation.hi_max_score, gene_info.hi_score)

    def _analyze_partial_overlap_loss(self, gene_info: GeneDosageInfo,
                                       annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """分析Loss部分重叠情况 - Section 2C/2D"""
        gene = gene_info.gene
        cnv = annotation.cnv

        # 获取该基因重叠的外显子
        gene_exons = [e for e in annotation.overlapping_exons if e.gene == gene]
        if not gene_exons:
            # 无外显子数据，使用默认评分
            evidence.section2_score = max(evidence.section2_score, 0.45)
            evidence.evidence_2B = 1
            evidence.evidence_details.append(
                f"Section 2B: Partial overlap with HI gene {gene} (no exon data, score=0.45)"
            )
            return

        # 检查是否涉及第一外显子(5'端)
        affects_first_exon = gene in annotation.affected_first_exons
        # 检查是否涉及最后外显子(3'端)
        affects_last_exon = gene in annotation.affected_last_exons

        # 统计重叠的外显子数量
        total_exons = self.db.get_gene_exon_count(gene)
        overlapped_count = len(gene_exons)

        if affects_first_exon:
            # Section 2C: 5'端部分重叠
            evidence.evidence_2C = 1
            score = 0.90  # 涉及编码序列
            evidence.section2_score = max(evidence.section2_score, score)
            evidence.evidence_details.append(
                f"Section 2C-1: Partial 5' overlap with {gene}, "
                f"affects first exon (exons overlapped: {overlapped_count}/{total_exons}, score={score:.2f})"
            )

        elif affects_last_exon:
            # Section 2D: 3'端部分重叠
            evidence.evidence_2D = 1
            # 如果只涉及最后外显子: 0.45-0.90; 多个外显子: 0.45-1.00
            if overlapped_count > 1:
                score = 0.90
                evidence.evidence_details.append(
                    f"Section 2D-4: Partial 3' overlap with {gene}, "
                    f"multiple exons (score={score:.2f})"
                )
            else:
                score = 0.45  # 仅最后外显子，无其他致病报道
                evidence.evidence_details.append(
                    f"Section 2D-3: Partial 3' overlap with {gene}, "
                    f"last exon only (score={score:.2f}) - MANUAL: check for other pathogenic variants"
                )
                evidence.manual_decision_needed = 1

            evidence.section2_score = max(evidence.section2_score, score)

        else:
            # 中间外显子或无法确定
            evidence.section2_score = max(evidence.section2_score, 0.45)
            evidence.evidence_2B = 1
            evidence.evidence_details.append(
                f"Section 2B: Partial overlap with {gene}, "
                f"internal exons (score=0.45)"
            )

    def _score_section2_gain(self, annotation: CNVAnnotation, evidence: ClinGenEvidence,
                             pathogenic_regions: List, overlapping_genes: List):
        """Gain (Duplication) Section 2评分 - 使用外显子信息"""
        cnv = annotation.cnv

        # 检查完全重叠TS score=3的基因
        for gene_info in overlapping_genes:
            if gene_info.tr_score == 3:
                if cnv.start <= gene_info.start and cnv.end >= gene_info.end:
                    evidence.section2_score = SECTION2_GAIN_SCORES["2A_full_overlap_TS"]
                    evidence.evidence_2A = 1
                    evidence.evidence_details.append(f"Section 2A: Full overlap with established TS gene {gene_info.gene} (score=1.00)")
                    annotation.tr_max_score = 3
                    return
                else:
                    # 部分重叠 - 需要分析外显子情况
                    self._analyze_partial_overlap_gain(gene_info, annotation, evidence)
                    annotation.tr_max_score = max(annotation.tr_max_score, gene_info.tr_score)

        # 检查致病区域重叠
        for region in pathogenic_regions:
            if region.tr_score == 3:
                overlap_pct = self._calculate_overlap_pct(cnv, region)
                if overlap_pct >= 0.80:
                    evidence.section2_score = SECTION2_GAIN_SCORES["2A_full_overlap_TS"]
                    evidence.evidence_2A = 1
                    evidence.evidence_details.append(f"Section 2A: High overlap with TS region {region.region_name} (score=1.00)")
                    annotation.tr_max_score = 3
                    return

        # 2I: 检查两断点是否都在HI基因内 (基因被中断)
        for gene_info in overlapping_genes:
            if gene_info.hi_score >= 2:  # HI score >= 2 (established or likely)
                # CNV是否完全包含在基因内 (两断点都在基因内)
                if gene_info.start <= cnv.start and cnv.end <= gene_info.end:
                    evidence.section2_score = max(evidence.section2_score, 0.90)
                    evidence.evidence_2E = 1  # 两断点在同一基因内(PVS1-like)
                    evidence.evidence_details.append(
                        f"Section 2I: Both breakpoints within HI gene {gene_info.gene} "
                        f"(HI={gene_info.hi_score}, gene disruption possible, score=0.90)"
                    )
                    annotation.hi_max_score = max(annotation.hi_max_score, gene_info.hi_score)

        # 2K: 一断点在HI基因内 + 表型一致 (需要人工决策)
        for gene_info in overlapping_genes:
            if gene_info.hi_score >= 2:
                # 一端在基因内，另一端在基因外
                if (gene_info.start <= cnv.start <= gene_info.end and cnv.end > gene_info.end) or \
                   (gene_info.start <= cnv.end <= gene_info.end and cnv.start < gene_info.start):
                    evidence.section2_score = max(evidence.section2_score, 0.45)
                    evidence.evidence_2K = 1  # 需要添加这个字段
                    evidence.evidence_details.append(
                        f"Section 2K: One breakpoint in HI gene {gene_info.gene} "
                        f"(HI={gene_info.hi_score}, score=0.45) - MANUAL: confirm phenotype consistency"
                    )
                    evidence.manual_decision_needed = 1
                    annotation.hi_max_score = max(annotation.hi_max_score, gene_info.hi_score)

        # TS score=2 的基因
        for gene_info in overlapping_genes:
            if gene_info.tr_score == 2:
                evidence.section2_score = max(evidence.section2_score, 0.15)
                evidence.evidence_2H = 1
                evidence.evidence_details.append(f"Section 2H: Likely TS gene {gene_info.gene} (score=0.15)")
                annotation.tr_max_score = max(annotation.tr_max_score, gene_info.tr_score)

    def _analyze_partial_overlap_gain(self, gene_info: GeneDosageInfo,
                                       annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """分析Gain部分重叠情况"""
        gene = gene_info.gene
        cnv = annotation.cnv

        # 获取该基因重叠的外显子
        gene_exons = [e for e in annotation.overlapping_exons if e.gene == gene]
        if not gene_exons:
            # 无外显子数据，使用默认评分
            evidence.section2_score = max(evidence.section2_score, 0.45)
            evidence.evidence_2B = 1
            evidence.evidence_details.append(
                f"Section 2B: Partial overlap with TS gene {gene} (no exon data, score=0.45)"
            )
            return

        # 统计重叠的外显子数量
        total_exons = self.db.get_gene_exon_count(gene)
        overlapped_count = len(gene_exons)

        # 检查是否涉及5'或3'端外显子
        affects_first_exon = gene in annotation.affected_first_exons
        affects_last_exon = gene in annotation.affected_last_exons

        if affects_first_exon or affects_last_exon:
            # 5'或3'端部分重叠
            evidence.section2_score = max(evidence.section2_score, 0.90)
            evidence.evidence_2C = 1
            evidence.exon_overlap_details.append(
                f"{gene}: partial overlap affects "
                f"{overlapped_count}/{total_exons} exons including 5'/3' end"
            )
            evidence.evidence_details.append(
                f"Section 2B: Partial overlap with TS gene {gene}, "
                f"affects terminal exon (exons overlapped: {overlapped_count}/{total_exons}, score=0.90)"
            )
        else:
            # 中间外显子
            evidence.section2_score = max(evidence.section2_score, 0.45)
            evidence.evidence_2B = 1
            evidence.evidence_details.append(
                f"Section 2B: Partial overlap with TS gene {gene}, "
                f"internal exons (score=0.45)"
            )

    def _score_section3(self, annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """Section 3: 基因数量评分"""
        cnv = annotation.cnv
        gene_count = annotation.gene_count

        if cnv.status == "Deletion":
            score_table = SECTION3_SCORES["loss"]
            thresholds = [50, 25, 10, 5, 2, 1, 0]
        else:
            score_table = SECTION3_SCORES["gain"]
            thresholds = [50, 35, 0]

        score = 0.0
        for threshold in thresholds:
            if gene_count >= threshold:
                score = score_table[threshold]
                break

        evidence.section3_score = score
        evidence.evidence_3_gene_count = 1 if score > 0 else 0
        evidence.evidence_details.append(f"Section 3: {gene_count} genes (score={score:.2f})")

    def _score_section4(self, annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """Section 4: 频率/数据库证据"""
        cnv = annotation.cnv

        # 查询频率
        freq_hits = self.db.query_frequency_by_interval(cnv.chrom, cnv.start, cnv.end, cnv.status)
        max_freq = 0.0

        for hit in freq_hits:
            if hit.get('source') == 'gnomAD':
                freq = hit.get('freq', 0)
                max_freq = max(max_freq, freq)
                annotation.frequency_info.append(FrequencyInfo(
                    source='gnomAD',
                    deletion_freq=freq if hit.get('svtype') == 'DEL' else 0,
                    duplication_freq=freq if hit.get('svtype') == 'DUP' else 0,
                    sample_size=0
                ))
            elif hit.get('source') == 'DECIPHER':
                del_freq = hit.get('del_freq', 0)
                dup_freq = hit.get('dup_freq', 0)
                if cnv.status == "Deletion":
                    max_freq = max(max_freq, del_freq)
                else:
                    max_freq = max(max_freq, dup_freq)
                annotation.frequency_info.append(FrequencyInfo(
                    source='DECIPHER',
                    deletion_freq=del_freq,
                    duplication_freq=dup_freq,
                    sample_size=hit.get('sample_size', 0)
                ))

        annotation.max_freq = max_freq

        # Section 4O: 频率证据
        if max_freq >= FREQUENCY_THRESHOLDS["benign"]:  # >= 1%
            evidence.section4_score = SECTION4_SCORES["4O_overlap_with_common_population"]
            evidence.evidence_4O = 1
            evidence.evidence_details.append(f"Section 4O: Frequency {max_freq:.2%} >= 1% (score=-1.00)")
        elif max_freq >= FREQUENCY_THRESHOLDS["likely_benign"]:  # >= 0.1%
            evidence.section4_score = -0.45  # 中等频率，部分良性证据
            evidence.evidence_4O = 1
            evidence.evidence_details.append(f"Section 4O-like: Frequency {max_freq:.2%} >= 0.1% (score=-0.45)")
        else:
            evidence.evidence_details.append(f"Section 4: Frequency {max_freq:.4f} < 0.1% (no benign frequency evidence)")

        # 标记不可用证据
        evidence.evidence_4A = 0
        evidence.evidence_details.append("Section 4A: No de novo data available")
        evidence.evidence_4L = 0
        evidence.evidence_details.append("Section 4L: No case-control data available")

    def _calculate_overlap_pct(self, cnv: CNVRecord, region) -> float:
        """计算重叠百分比"""
        overlap_start = max(cnv.start, region.start)
        overlap_end = min(cnv.end, region.end)
        overlap_len = overlap_end - overlap_start
        region_len = region.end - region.start
        return overlap_len / region_len if region_len > 0 else 0

    def _classify(self, annotation: CNVAnnotation, evidence: ClinGenEvidence):
        """最终分类"""
        total = evidence.total_score

        if total <= CLASSIFICATION_THRESHOLDS["Benign"]:
            annotation.overall_classification = "Benign"
            annotation.classification_reason = f"Total score {total:.2f} <= -0.99"
        elif total <= CLASSIFICATION_THRESHOLDS["VUS_upper"]:
            annotation.overall_classification = "Uncertain Significance"
            annotation.classification_reason = f"Total score {total:.2f} in VUS range (-0.98 to 0.99)"
        elif total >= CLASSIFICATION_THRESHOLDS["Pathogenic"]:
            annotation.overall_classification = "Pathogenic"
            annotation.classification_reason = f"Total score {total:.2f} >= 1.00"
        else:
            annotation.overall_classification = "Likely Pathogenic"
            annotation.classification_reason = f"Total score {total:.2f}"


# ============================================================================
# CNV注释器
# ============================================================================

class CNVAnnotator:
    """CNV注释器"""

    def __init__(self, db_loader: DatabaseLoader):
        self.db = db_loader
        self.scorer = ClinGenScorer(db_loader)
        self.iscn_formatter = ISCNFormatter(db_loader.genome_build)
        self.logger = logging.getLogger(__name__)

    def annotate_cnv(self, cnv: CNVRecord) -> CNVAnnotation:
        """注释单个CNV"""
        annotation = CNVAnnotation(cnv=cnv)

        if cnv.status == "Normal":
            annotation.overall_classification = "Benign"
            annotation.classification_reason = "No CNV detected"
            annotation.iscn = self.iscn_formatter.format_iscn(
                cnv.chrom, cnv.start, cnv.end, cnv.status, "", ""
            )
            return annotation

        # 执行ClinGen评分
        annotation.clingen_evidence = self.scorer.score_cnv(annotation)

        # 补充GenCC信息
        self._annotate_gencc(annotation)

        # 生成ISCN标准命名
        start_cytoband = self.db.get_cytoband(cnv.chrom, cnv.start)
        end_cytoband = self.db.get_cytoband(cnv.chrom, cnv.end)
        annotation.iscn = self.iscn_formatter.format_iscn(
            cnv.chrom, cnv.start, cnv.end, cnv.status, start_cytoband, end_cytoband
        )

        return annotation

    def _annotate_gencc(self, annotation: CNVAnnotation):
        """补充GenCC基因-疾病信息"""
        supportive_genes = set()
        for gene in annotation.overlapping_genes:
            if gene in self.db.gencc_data:
                for gencc_info in self.db.gencc_data[gene]:
                    annotation.gencc_info.append(gencc_info)
                    # AD遗传模式的基因对CNV更有意义
                    if gencc_info.moi in ["Autosomal dominant", "X-linked"]:
                        supportive_genes.add(gene)
        annotation.gencc_supportive_genes = sorted(supportive_genes)


# ============================================================================
# 输入解析器
# ============================================================================

class InputParser:
    """解析输入CNV文件"""

    def __init__(self, genome_build: str = "GRCh38"):
        self.genome_build = genome_build
        self.logger = logging.getLogger(__name__)

    def parse_file(self, input_file: str) -> List[CNVRecord]:
        """解析输入文件"""
        path = Path(input_file)

        if path.suffix == '.cnr':
            cnvs = self._parse_cnr(input_file)
        elif path.suffix in ['.bed', '.txt', '.tsv']:
            cnvs = self._parse_bed(input_file)
        elif path.suffix == '.vcf' or str(path).endswith('.vcf.gz'):
            cnvs = self._parse_vcf(input_file)
        else:
            cnvs = self._parse_generic(input_file)

        self.logger.info(f"Parsed {len(cnvs)} CNV records")
        return cnvs

    def _parse_cnr(self, input_file: str) -> List[CNVRecord]:
        """解析CNVkit .cnr格式"""
        cnvs = []
        df = pd.read_csv(input_file, sep='\t')

        for _, row in df.iterrows():
            chrom = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            log2 = float(row['log2']) if 'log2' in row else 0

            if log2 < -0.6:
                status = "Deletion"
            elif log2 > 0.4:
                status = "Amplification"
            else:
                status = "Normal"

            cnvs.append(CNVRecord(chrom=chrom, start=start, end=end, status=status))

        return cnvs

    def _parse_bed(self, input_file: str) -> List[CNVRecord]:
        """解析BED格式(4列)"""
        cnvs = []

        with open(input_file) as f:
            for line in f:
                if line.startswith('#') or line.startswith('track') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                if parts[1] == 'start' or parts[2] == 'end' or parts[3] == 'status':
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                status = parts[3].strip()

                cnvs.append(CNVRecord(chrom=chrom, start=start, end=end, status=status))

        return cnvs

    def _parse_vcf(self, input_file: str) -> List[CNVRecord]:
        """解析VCF格式"""
        cnvs = []
        opener = gzip.open if input_file.endswith('.gz') else open
        mode = 'rt' if input_file.endswith('.gz') else 'r'

        with opener(input_file, mode) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue

                chrom = parts[0]
                pos = int(parts[1])

                info_dict = {}
                for item in parts[7].split(';'):
                    if '=' in item:
                        key, val = item.split('=')
                        info_dict[key] = val

                end = int(info_dict.get('END', pos))
                svtype = info_dict.get('SVTYPE', '')

                if svtype == 'DEL':
                    status = "Deletion"
                elif svtype == 'DUP':
                    status = "Amplification"
                else:
                    status = "Normal"

                cnvs.append(CNVRecord(chrom=chrom, start=pos, end=end, status=status))

        return cnvs

    def _parse_generic(self, input_file: str) -> List[CNVRecord]:
        """通用格式解析"""
        df = pd.read_csv(input_file, sep='\t')

        chrom_col = self._find_column(df, ['chrom', 'chromosome', 'chr'])
        start_col = self._find_column(df, ['start', 'begin'])
        end_col = self._find_column(df, ['end', 'stop'])
        status_col = self._find_column(df, ['status', 'type', 'cnv_type'])

        if not all([chrom_col, start_col, end_col]):
            raise ValueError("Cannot identify required columns")

        cnvs = []
        for _, row in df.iterrows():
            status = row[status_col] if status_col else "Uncertain"
            cnvs.append(CNVRecord(
                chrom=str(row[chrom_col]),
                start=int(row[start_col]),
                end=int(row[end_col]),
                status=str(status)
            ))

        return cnvs

    def _find_column(self, df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
        """查找匹配的列名"""
        for col in df.columns:
            col_lower = col.lower()
            for candidate in candidates:
                if candidate in col_lower:
                    return col
        return None


# ============================================================================
# 输出格式化
# ============================================================================

class OutputFormatter:
    """格式化输出结果"""

    def __init__(self, output_format: str = "tsv", db_loader=None):
        self.output_format = output_format
        self.db = db_loader

    def format_results(self, annotations: List[CNVAnnotation]) -> str:
        if self.output_format == "tsv":
            return self._format_tsv(annotations)
        elif self.output_format == "json":
            return self._format_json(annotations)
        else:
            return self._format_tsv(annotations)

    def _format_tsv(self, annotations: List[CNVAnnotation]) -> str:
        """TSV格式输出"""
        headers = [
            "#Chromosome", "Start", "End", "Size", "Status", "ISCN",
            "Gene_Count", "HI_Max", "TR_Max", "Max_Frequency",
            # ClinGen评分
            "Section1", "Section2", "Section3", "Section4", "Section5", "Total_Score",
            # 证据标记
            "Evidence_1A", "Evidence_1B", "Evidence_2A", "Evidence_2B",
            "Evidence_2C", "Evidence_2D", "Evidence_2E", "Evidence_2F", "Evidence_2H", "Evidence_2K",
            "Evidence_3", "Evidence_4O", "Evidence_4A", "Evidence_4L", "Evidence_5",
            # 其他信息
            "Dosage_Genes", "Pathogenic_Regions", "Benign_Regions_Overlap",
            "GenCC_AD_Genes", "Classification", "Reason", "Evidence_Details"
        ]

        lines = ['\t'.join(headers)]

        for ann in annotations:
            cnv = ann.cnv
            ev = ann.clingen_evidence

            # 格式化基因和区域
            ds_genes = ';'.join([f"{g.gene}(HI={g.hi_score},TR={g.tr_score})"
                                 for g in ann.dosage_sensitive_genes[:10]])
            path_regions = ';'.join([r.region_name for r in ann.overlapping_regions[:5]])

            # 检查良性区域重叠
            benign_overlap = ""
            if self.db:
                benign_regions = self.db.query_benign_regions_by_interval(cnv.chrom, cnv.start, cnv.end)
                benign_overlap = ';'.join([r.region_name for r in benign_regions[:3]])

            row = [
                cnv.chrom, cnv.start, cnv.end, cnv.size, cnv.status, ann.iscn,
                ann.gene_count, ann.hi_max_score, ann.tr_max_score,
                f"{ann.max_freq:.6f}",
                # 评分
                f"{ev.section1_score:.2f}", f"{ev.section2_score:.2f}",
                f"{ev.section3_score:.2f}", f"{ev.section4_score:.2f}",
                f"{ev.section5_score:.2f}", f"{ev.total_score:.2f}",
                # 证据标记
                ev.evidence_1A, ev.evidence_1B, ev.evidence_2A, ev.evidence_2B,
                ev.evidence_2C, ev.evidence_2D, ev.evidence_2E, ev.evidence_2F,
                ev.evidence_2H, ev.evidence_2K, ev.evidence_3_gene_count,
                ev.evidence_4O, ev.evidence_4A, ev.evidence_4L,
                ev.evidence_5,
                # 其他
                ds_genes, path_regions, benign_overlap,
                ','.join(ann.gencc_supportive_genes[:10]),
                ann.overall_classification,
                ann.classification_reason,
                '; '.join(ev.evidence_details[:20])
            ]

            lines.append('\t'.join(str(x) for x in row))

        return '\n'.join(lines)

    def _format_json(self, annotations: List[CNVAnnotation]) -> str:
        """JSON格式输出"""
        results = []
        for ann in annotations:
            results.append({
                "cnv": {
                    "chromosome": ann.cnv.chrom,
                    "start": ann.cnv.start,
                    "end": ann.cnv.end,
                    "size": ann.cnv.size,
                    "status": ann.cnv.status,
                    "iscn": ann.iscn
                },
                "annotation": {
                    "gene_count": ann.gene_count,
                    "dosage_genes": [g.gene for g in ann.dosage_sensitive_genes],
                    "hi_max": ann.hi_max_score,
                    "tr_max": ann.tr_max_score,
                    "frequency": ann.max_freq,
                    "clingen_score": {
                        "section1": ann.clingen_evidence.section1_score,
                        "section2": ann.clingen_evidence.section2_score,
                        "section3": ann.clingen_evidence.section3_score,
                        "section4": ann.clingen_evidence.section4_score,
                        "section5": ann.clingen_evidence.section5_score,
                        "total": ann.clingen_evidence.total_score
                    },
                    "evidence_used": {
                        "1A_protein_coding": ann.clingen_evidence.evidence_1A,
                        "1B_no_protein_coding": ann.clingen_evidence.evidence_1B,
                        "2A_full_overlap_pathogenic": ann.clingen_evidence.evidence_2A,
                        "2B_partial_overlap": ann.clingen_evidence.evidence_2B,
                        "2C_5prime_partial": ann.clingen_evidence.evidence_2C,
                        "2D_3prime_partial": ann.clingen_evidence.evidence_2D,
                        "2E_both_breakpoints_in_gene": ann.clingen_evidence.evidence_2E,
                        "2F_benign_region": ann.clingen_evidence.evidence_2F,
                        "2H_HI_predictor": ann.clingen_evidence.evidence_2H,
                        "2K_one_breakpoint_in_HI": ann.clingen_evidence.evidence_2K,
                        "3_gene_count": ann.clingen_evidence.evidence_3_gene_count,
                        "4O_high_frequency": ann.clingen_evidence.evidence_4O,
                        "4A_de_novo": ann.clingen_evidence.evidence_4A,
                        "4L_case_control": ann.clingen_evidence.evidence_4L,
                        "5_inheritance": ann.clingen_evidence.evidence_5
                    },
                    "exon_overlap_details": ann.clingen_evidence.exon_overlap_details,
                    "manual_decision_needed": ann.clingen_evidence.manual_decision_needed,
                    "evidence_details": ann.clingen_evidence.evidence_details,
                    "classification": ann.overall_classification,
                    "reason": ann.classification_reason
                }
            })
        return json.dumps(results, indent=2)

    def write_output(self, annotations: List[CNVAnnotation], output_file: str):
        content = self.format_results(annotations)
        with open(output_file, 'w') as f:
            f.write(content)


# ============================================================================
# 主程序
# ============================================================================

def setup_logging(log_level: str = "INFO"):
    level = getattr(logging, log_level.upper())
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def main():
    parser = argparse.ArgumentParser(
        description='CNVAnno - ClinGen Official CNV Pathogenicity Scorer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
ClinGen Official Scoring System:
  Section 1: Genomic content (0 or -0.60)
  Section 2: Gene/region overlap (0 to 1.00, or -1.00 for benign)
  Section 3: Gene count (0 to 0.90)
  Section 4: Frequency/database evidence (-1.00 for high freq)
  Section 5: Inheritance (data unavailable = 0)

Classification:
  Total <= -0.99: Benign
  Total -0.98 ~ 0.99: VUS
  Total >= 1.00: Pathogenic

Usage:
  cnvanno input.bed -o output.tsv
  cnvanno input.bed -g GRCh37 -o output.tsv
        """
    )

    parser.add_argument('input', help='Input CNV file (4 columns: chrom, start, end, status)')
    parser.add_argument('-d', '--data-dir', default='data', help='Database directory')
    parser.add_argument('-g', '--genome-build', default='GRCh38', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-f', '--format', default='tsv', choices=['tsv', 'json'])
    parser.add_argument('-l', '--log-level', default='INFO')

    args = parser.parse_args()
    setup_logging(args.log_level)
    logger = logging.getLogger('CNVAnno')

    # 检查数据目录
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        sys.exit(1)

    # 加载数据库
    logger.info("Loading ClinGen databases...")
    db_loader = DatabaseLoader(args.data_dir, args.genome_build)
    db_loader.load_all()

    # 解析输入
    logger.info(f"Parsing input: {args.input}")
    input_parser = InputParser(args.genome_build)
    cnvs = input_parser.parse_file(args.input)

    # 注释CNV
    logger.info("Annotating CNVs with ClinGen scoring...")
    annotator = CNVAnnotator(db_loader)
    annotations = [annotator.annotate_cnv(cnv) for cnv in cnvs]

    # 输出结果
    logger.info(f"Writing output: {args.output}")
    formatter = OutputFormatter(args.format, db_loader)
    formatter.write_output(annotations, args.output)

    # 统计
    classifications = defaultdict(int)
    for ann in annotations:
        classifications[ann.overall_classification] += 1

    logger.info("Classification summary:")
    for cls, count in sorted(classifications.items()):
        logger.info(f"  {cls}: {count}")

    logger.info("Done!")


if __name__ == '__main__':
    main()