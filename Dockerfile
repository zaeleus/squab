FROM rust:1.41.0-buster as builder

COPY Cargo.lock Cargo.toml /tmp/noodles-squab/
COPY src/ /tmp/noodles-squab/src/

RUN cargo build \
      --release \
      --manifest-path /tmp/noodles-squab/Cargo.toml

FROM debian:buster

COPY --from=builder \
    /tmp/noodles-squab/target/release/noodles-squab \
    /opt/noodles-squab/bin/

ENTRYPOINT ["/opt/noodles-squab/bin/noodles-squab"]
