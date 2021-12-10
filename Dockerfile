FROM rust:1.57.0-bullseye as builder

WORKDIR /tmp/noodles-squab/

COPY .git/ .git/
COPY Cargo.lock Cargo.toml ./
COPY src/ src/

RUN cargo build --release

FROM debian:bullseye

COPY --from=builder \
    /tmp/noodles-squab/target/release/noodles-squab \
    /opt/noodles-squab/bin/

ENV PATH=/opt/noodles-squab/bin:$PATH

ENTRYPOINT ["/opt/noodles-squab/bin/noodles-squab"]
