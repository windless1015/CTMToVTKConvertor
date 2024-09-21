cmake_minimum_required(VERSION 3.5)

project(CTMToVTKConvertor VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#include directory
include_directories(common)

#qt setting
find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets)


#vtk library setting
find_package(VTK REQUIRED)

#openctm library setting
include_directories(thirdParty) #搜索头文件路径
link_directories(thirdParty) #添加需要链接的库文件路径




set(PROJECT_SOURCES
    main.cpp
	CTMIO.h
	CTMIO.cpp
)

add_executable(CTMToVTKConvertor ${PROJECT_SOURCES}) #向项目添加可执行文件构建目标
target_link_libraries(CTMToVTKConvertor PRIVATE 
Qt${QT_VERSION_MAJOR}::Widgets
${VTK_LIBRARIES}
openctm.lib
)

set_target_properties(CTMToVTKConvertor PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER my.example.com
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

