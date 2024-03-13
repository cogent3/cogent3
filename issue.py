from cogent3 import open_data_store
ds = open_data_store("loci.zip", suffix="phy", mode="r")
ds.members[:5]
print(ds.members[:5])