FROM python:3.11

WORKDIR /usr/src/ct_dataset

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

ENTRYPOINT ["python"]