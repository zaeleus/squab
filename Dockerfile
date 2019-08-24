FROM rust:1.37.0-buster as builder

COPY Cargo.lock Cargo.toml /tmp/noodles-count-features/
COPY src/ /tmp/noodles-count-features/src/

RUN cargo build \
      --release \
      --manifest-path /tmp/noodles-count-features/Cargo.toml

FROM debian:buster

COPY --from=builder \
    /tmp/noodles-count-features/target/release/noodles-count-features \
    /opt/noodles-count-features/bin/

ENTRYPOINT ["/opt/noodles-count-features/bin/noodles-count-features"]
