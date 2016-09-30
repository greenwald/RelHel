include_directories(${RelHel_SOURCE_DIR}/include)

set(PROGRAMS
  test
#  CalcAmpl
  )

foreach(program ${PROGRAMS})
	add_executable(${program} ${program}.cxx)
	target_link_libraries(${program} RelHel)
endforeach(program)
