# KMEA
> [!NOTE]
> This is an academic project, the instructions can be found under the doc directory (or directly [here](./doc/INFOH417%20DBSA%20–%20Project%2024-25.pdf))


KMEA[^1] is a PostgreSQL extension that supports various DNA data types, along with some operators.
## Supported data types
- DNA sequences
- Kmers
- Qkmers

## Available functions
- Length
- Startswith
- Equals
- Qkmer contains Kmer
- Generate Kmers

## Additional features
- Hash function for kmer counting support
- SP-GiST index for kmers


[^1]: Kmer Extension for Analysis
---

# How to setup
> [!IMPORTANT]
> Make sure you have the following installed:
>   - PostgreSQL 16
>   - PostgreSQL server dev 16

Make the extension with the following commands
```shell
make
sudo make install
```
---
# Testing features
You can either create the extension and test by yourself
```shell
createdb kmea
psql kmea
```
```sql
create extension kmea;
```
Or launch the [kmea test db file](kmea_test_db.sql)
```shell
psql -d postgres -f kmea_test_db.sql
```
