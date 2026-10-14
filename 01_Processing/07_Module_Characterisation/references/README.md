# Module characterisation references

`STRING/` holds the human STRING v12.0 files read by `../a_script/module_characterisation.R`:

| File | Purpose |
|---|---|
| `9606.protein.aliases.v12.0.txt.gz` | Map protein accessions and alternative names to STRING identifiers |
| `9606.protein.info.v12.0.txt.gz` | Protein names and annotations |
| `9606.protein.links.v12.0.min700.txt.gz` | Connections with combined score at least 700 |

Source: <https://string-db.org/cgi/download?species_text=Homo+sapiens>, version 12.0.
The script uses these local files for identifier mapping and network analysis. Retain the
versioned files together; no new download step is introduced by this layout.
