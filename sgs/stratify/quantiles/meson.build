if WIN
  #create breaks module
  py.extension_module(
    'quantiles',
    'quantiles.cpp',
    link_with: [
      raster_lib],
    install: true,
  )
endif
