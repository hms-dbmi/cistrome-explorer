After running the notebook, to ingest the tileset using higlass server's `manage.py` script:

```sh
python manage.py ingest_tileset \
    --uid 'Cistrome_DNase_1kb_average' \
    --filename ./Cistrome_DNase_1kb_average_chr1_to_chr22.multires.mv5 \
    --filetype multivec \
    --coordSystem hg38

```