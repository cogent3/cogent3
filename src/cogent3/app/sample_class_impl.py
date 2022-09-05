import pluggy
from cogent3.app.composable import define_app

hookimpl = pluggy.HookimplMarker(project_name="datastore")

@define_app
class sample_app:
    
    def __init__(self):
        pass

    @hookimpl
    def main(self, a: int) -> int:
        return a * 7
