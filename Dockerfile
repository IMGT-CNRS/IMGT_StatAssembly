FROM rust:latest

# Install zip utility (needed for creating the archive)
RUN apt-get update && apt-get install -y zip libclang-dev libssl-dev && rm -rf /var/lib/apt/lists/*


WORKDIR /app
RUN git clone https://github.com/IMGT-CNRS/IMGT_StatAssembly.git . && \
    rm -rf .git && \
    cargo build --release && \
    mkdir -p /output && \
    cp target/release/IMGT_StatAssembly /output/ && \
    cp README* /output/ 2>/dev/null || true && \
    cp LICENSE* /output/ 2>/dev/null || true && \
    cd /output && \
    zip IMGT_StatAssembly.zip IMGT_StatAssembly README* LICENSE* 2>/dev/null || true
# Create the zip archive

# Set the default command to print instructions (or override at runtime)
CMD ["sh", "-c", "cp /output/IMGT_StatAssembly.zip /host/ && echo 'Zip file copied to host output folder.'"]
