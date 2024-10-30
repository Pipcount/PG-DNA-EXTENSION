# KMEA extension

## How to setup
1. Make sure you have the following installed:
   - PostgreSQL
   - PostgreSQL server
2. Create a database in PostgreSQL
```shell
createdb kmea
```
3. Make the extension:
```shell
make
sudo make install
```

4. Create the extension:
```shell
psql kmea
```
```sql
create extension kmea;
```