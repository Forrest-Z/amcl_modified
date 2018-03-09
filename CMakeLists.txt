cmake_minimum_required(VERSION 2.8.3)
project(amcl_modified)

find_package(catkin REQUIRED
        COMPONENTS
            message_filters
            rosbag
            roscpp
            std_srvs
            tf
            dynamic_reconfigure
            nav_msgs
        )

find_package(Boost REQUIRED)

# dynamic reconfigure
generate_dynamic_reconfigure_options(
    cfg/AMCL.cfg
)

catkin_package(
    CATKIN_DEPENDS
        rosbag
        roscpp
        dynamic_reconfigure
        tf
  CATKIN_DEPENDS nav_msgs std_srvs
  INCLUDE_DIRS include
  LIBRARIES amcl_modified_sensors amcl_modified_map amcl_modified_pf
)

include_directories(include)
include_directories(${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_library(amcl_modified_pf
                    src/amcl/pf/pf.c
                    src/amcl/pf/pf_kdtree.c
                    src/amcl/pf/pf_pdf.c
                    src/amcl/pf/pf_vector.c
                    src/amcl/pf/eig3.c
                    src/amcl/pf/pf_draw.c)

add_library(amcl_modified_map
                    src/amcl/map/map.c
                    src/amcl/map/map_cspace.cpp
                    src/amcl/map/map_range.c
                    src/amcl/map/map_store.c
                    src/amcl/map/map_draw.c)

add_library(amcl_modified_sensors
                    src/amcl/sensors/amcl_sensor.cpp
                    src/amcl/sensors/amcl_odom.cpp
                    src/amcl/sensors/amcl_laser.cpp)
target_link_libraries(amcl_modified_sensors amcl_modified_map amcl_modified_pf)

install( TARGETS
    amcl_modified_sensors amcl_modified_map amcl_modified_pf
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/amcl/
  DESTINATION "include/amcl"
)

install(DIRECTORY examples/
    DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/examples
)
