### docker build
docker build -t ct_dataset .

### docker run
Will download ChEMBL every time unless SQLite path is given
docker run -it --rm -v "$PWD":/usr/src/ct_dataset ct_dataset main.py --help
docker run -it --rm -v "$PWD":/usr/src/ct_dataset ct_dataset main.py -v 32 -o OUTPUT_PATH

