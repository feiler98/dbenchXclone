FROM python:3.9

RUN mkdir -p /scratch/tmp/feiler/dbenchXclone
WORKDIR /scratch/tmp/feiler/dbenchXclone
COPY . .

RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/scratch/tmp/feiler/dbenchXclonerun_xclone.py"]