FROM python:buster


# Install system dependencies
RUN git clone https://github.com/USDA-ARS-GBRU/micrompn.git && \
     cd micrompn && \
     pip install --no-cache-dir .


# Set the default command to run itsxpress
CMD ["micrompn"]