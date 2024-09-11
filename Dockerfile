FROM rust:1.81.0 AS builder

COPY . /psdm

WORKDIR /psdm

ARG TARGET="x86_64-unknown-linux-musl"

RUN apt update \
    && apt install -y musl-tools \
    && rustup target add "$TARGET" \
    && cargo build --release --target "$TARGET" \
    && strip target/${TARGET}/release/psdm


FROM bash:5.0

ARG TARGET="x86_64-unknown-linux-musl"
COPY --from=builder /psdm/target/${TARGET}/release/psdm /bin/

RUN psdm --version

ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
