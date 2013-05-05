SET(vocabulary_SRCS 
  Vocabulary.cpp
) 

SET(vocabulary_HDRS 
  HierarchicalKMeansClustering.h
  HierarchicalKMeansClustering.hpp
  KMeansClustering.h
  KMeansClustering.hpp
  Vocabulary.h
) 

ADD_LIBRARY(vocabulary SHARED ${vocabulary_SRCS})
TARGET_LINK_LIBRARIES(vocabulary feature geometry)

install(TARGETS vocabulary
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib/${PROJECT_NAME}
    ARCHIVE DESTINATION lib/${PROJECT_NAME})

install(FILES ${vocabulary_HDRS} DESTINATION include/${PROJECT_NAME}/vocabulary)
