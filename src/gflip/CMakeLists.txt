SET(gflip_SRCS 
  gflip_engine.cpp
) 

SET(gflip_HDRS 
  gflip_engine.hpp
) 

ADD_LIBRARY(gflip SHARED ${gflip_SRCS})

install(TARGETS gflip
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib/${PROJECT_NAME}
    ARCHIVE DESTINATION lib/${PROJECT_NAME})

install(FILES ${gflip_HDRS} DESTINATION include/${PROJECT_NAME}/gflip)
