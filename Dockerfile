ARG RUST_VERSION=1.78.0
ARG DEBIAN_CODENAME=bookworm

FROM rust:${RUST_VERSION}-${DEBIAN_CODENAME} as builder

WORKDIR /tmp/squab/

COPY .git/ .git/
COPY Cargo.lock Cargo.toml ./
COPY src/ src/

RUN cargo build --release

FROM debian:${DEBIAN_CODENAME}

COPY --from=builder /tmp/squab/target/release/squab /opt/squab/bin/

ENV PATH=/opt/squab/bin:$PATH

ENTRYPOINT ["/opt/squab/bin/squab"]
