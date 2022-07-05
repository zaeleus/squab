FROM rust:1.62.0-bullseye as builder

WORKDIR /tmp/squab/

COPY .git/ .git/
COPY Cargo.lock Cargo.toml ./
COPY src/ src/

RUN cargo build --release

FROM debian:bullseye

COPY --from=builder /tmp/squab/target/release/squab /opt/squab/bin/

ENV PATH=/opt/squab/bin:$PATH

ENTRYPOINT ["/opt/squab/bin/squab"]
