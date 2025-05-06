# Use Ubuntu as base image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ARG http_proxy
ARG https_proxy

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    unzip \
    curl \
    gnupg \
    ca-certificates \
    build-essential \
    r-base \
    xvfb \
    libxi6 \
    libgconf-2-4 \
    libnss3 \
    libxss1 \
    libappindicator3-1 \
    fonts-liberation \
    libatk-bridge2.0-0 \
    libatk1.0-0 \
    libgtk-3-0 \
    libx11-xcb1 \
    libasound2 \
    libgbm1 \
    libvulkan1 \
    libnspr4 \
    xdg-utils \
    lsb-release \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages (e.g., ggplot2)
RUN Rscript -e "install.packages(c('ggplot2'), repos='http://cran.r-project.org')"

# Install Google Chrome (v136.0.7103.49)
RUN wget -O chrome-linux64.zip https://storage.googleapis.com/chrome-for-testing-public/136.0.7103.49/linux64/chrome-linux64.zip && \
    unzip chrome-linux64.zip && \
    mv chrome-linux64 /opt/google-chrome && \
    ln -s /opt/google-chrome/chrome /usr/bin/google-chrome && \
    rm chrome-linux64.zip

# Install ChromeDriver (v136.0.7103.49)
RUN wget -O chromedriver-linux64.zip https://storage.googleapis.com/chrome-for-testing-public/136.0.7103.49/linux64/chromedriver-linux64.zip && \
    unzip chromedriver-linux64.zip && \
    mv chromedriver-linux64/chromedriver /usr/local/bin/chromedriver && \
    chmod +x /usr/local/bin/chromedriver && \
    rm -rf chromedriver-linux64*

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    requests \
    biopython \
    pandas \
    selenium \
    google-generativeai \
    google-api-core    

# Set working directory and copy files
WORKDIR /home/eukprogs
COPY . /home/eukprogs
RUN mkdir -p /home/eukprogs/downloads
RUN mkdir -p /home/eukprogs/Results
RUN rm -rf /tmp/*

# Default command to run
CMD ["python3", "chiranth.py"]

