# KMEA extension
> [!NOTE]
> This is an academic project, the instructions can be found under the doc directory (or directly [here](./doc/INFOH417%20DBSA%20–%20Project%2024-25.pdf))


The KMEA extension is a PostgreSQL extension that supports various DNA data types, along with some operators.
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





---
# How to setup
> [!WARNING]
> Make sure you have the following installed:
>   - PostgreSQL 16
>   - PostgreSQL server dev 16
1. Create a database in PostgreSQL
```shell
createdb kmea
```
2. Make the extension:
```shell
make
sudo make install
```

3. Create the extension:
```shell
psql kmea
```
---
# Testing features
You can either create the extension and test by yourself
```sql
create extension kmea;
```
Or launch the [kmea test db file](kmea_test_db.sql)
```sql
\i kmea_test_db.sql
```