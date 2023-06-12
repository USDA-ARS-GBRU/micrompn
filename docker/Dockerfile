FROM python:buster
LABEL author="Adam Rivers"
LABEL description="MicroMPN: Software for automating most probable number (MPN) estimates from laboratory microplates"

# Install system dependencies
RUN git clone https://github.com/USDA-ARS-GBRU/micrompn.git && \
     cd micrompn && \
     pip install --no-cache-dir .


# Set the default command to run itsxpress
CMD ["micrompn"]