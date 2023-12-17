FROM python:3.11

WORKDIR /usr/src/ct_dataset

COPY pyproject.toml ./
RUN pip install .

ENTRYPOINT ["python"]