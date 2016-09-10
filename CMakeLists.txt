#generate coredll 


INCLUDE_DIRECTORIES(/usr/local/include  
										../ann_1.1_char
										)


ADD_LIBRARY(CoreDll
	main.cpp
	vim_imgs.cpp
	)