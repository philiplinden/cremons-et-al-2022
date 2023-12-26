FROM python:3.11

USER root

COPY pyproject.toml /opt
WORKDIR /opt

RUN pip install .
COPY . /opt

ENTRYPOINT ["python"]
