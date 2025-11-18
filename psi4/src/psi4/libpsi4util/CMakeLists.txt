list(APPEND sources
  PsiOutStream.cc
  combinations.cc
  exception.cc
  header_printer.cc
  memory_manager.cc
  process.cc
  stl_string.cc
  )
psi4_add_module(lib psi4util sources)
