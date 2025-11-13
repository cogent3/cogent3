from cogent3.app.comp_new import define_app


@define_app
class app_dummyclass_1:
    def __init__(self, a):
        self.a = a

    def main(self, val: int) -> int:
        return val


@define_app
class app_dummyclass_2:
    def __init__(self, b):
        self.b = b

    def main(self, val: int) -> int:
        return val * 2


@define_app
class app_dummyclass_3:
    def __init__(self, c):
        self.b = c

    def main(self, val: int) -> int:
        return val


aseqfunc1 = app_dummyclass_1(1)
aseqfunc2 = app_dummyclass_2(2)
aseqfunc3 = app_dummyclass_3(3)
comb = (aseqfunc2 + aseqfunc2) + aseqfunc2

print(comb(3))
