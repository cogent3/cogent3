
import pluggy

hookspec1 = pluggy.HookspecMarker(project_name="datastore")

class MyClassSpec:
    """A hook specification namespace."""

    @hookspec1
    def main(self, a:int):
        ...