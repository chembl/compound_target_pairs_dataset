# Generating / Updating Documentation

```
sphinx-apidoc -o ./source ../src
```

# Building Documentation
```
make clean
make html
make latex && cd build/latex && make
```
