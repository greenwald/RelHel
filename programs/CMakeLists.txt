include_directories(${RelHel_SOURCE_DIR}/include)

set(PROGRAMS CalcAmpl)
foreach(program ${PROGRAMS})
	add_executable(${program} EXCLUDE_FROM_ALL ${program}.cc)
	target_link_libraries(${program} RelHel)
endforeach(program)