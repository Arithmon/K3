# Reference environment for full reproduction.
#
# Every version below is derived from the libraries the chain actually
# imports, not from a wish list. The interval endpoints of the certificates
# depend on the arithmetic backend, so the backend is part of the claim.
FROM python:3.14.4-slim

ENV PYTHONHASHSEED=0 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    SOURCE_DATE_EPOCH=0

WORKDIR /k3
COPY requirements.lock ./
RUN pip install --no-cache-dir -r requirements.lock
COPY . .
RUN pip install --no-cache-dir -e .

CMD ["python", "verification/verify.py"]
