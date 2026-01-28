using TestItems

using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
