
# Recurse through the directory
file(GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

# Add to sources
set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
