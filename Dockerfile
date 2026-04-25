FROM python:3.11-slim

LABEL maintainer="CNVAnno"
LABEL description="ClinGen CNV Pathogenicity Annotation Tool"
LABEL version="1.0"

# Set working directory
WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY cnvanno.py .

# Create data directory (for mounting)
RUN mkdir -p /app/data

# Set entrypoint
ENTRYPOINT ["python", "cnvanno.py"]

# Default command shows help
CMD ["--help"]