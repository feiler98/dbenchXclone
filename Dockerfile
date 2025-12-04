FROM python:3.9

COPY . .

RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "run_xclone.py"]