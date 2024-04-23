Note: To avoid the mesmer module re-downloading the 100 MB model weights file
for every run of every test case, run nf-test with the `keras_cache_tmp`
profile. This will map /tmp/.keras on your physical filesystem into the
container to provide a durable cache. Use the `+profile` syntax to also retain
the default `docker` profile rather than replace it:

```shell
nf-test test --profile +keras_cache_tmp
```
