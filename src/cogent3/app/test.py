from cogent3.app.sample_class_spec import MyClassSpec as hookspec
from cogent3.app.sample_class_impl import sample_app as hookimpl
import pluggy

print("ttt")

def get_plugin_manager():
    pm = pluggy.PluginManager(project_name="datastore")
    pm.add_hookspecs(hookspec)
    pm.register(hookimpl())
    pm.load_setuptools_entrypoints("datastore")
    return pm


pm = get_plugin_manager()

#supports only keyword arguments
print( pm.hook.main(a=11) )

