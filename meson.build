project('purelib-and-platlib', 'c')

py = import('python').find_installation(pure: false)

py.extension_module(
    'test_module',
    'test_module.c',
    install: true,
)
