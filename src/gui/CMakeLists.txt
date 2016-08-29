#
# Qt user interface cmake values update to Qt5 (wip)
#
include_directories(${YAF_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if (EMBED_FONT_QT)
    add_definitions(-DYAFQT_EMBEDED_FONT)
endif (EMBED_FONT_QT)

add_definitions(-DBUILDING_QTPLUGIN)

set(YAF_QT_CC animworking.cc events.cc mywindow.cc qtoutput.cc qtprogress.cc renderwidget.cc worker.cc)
set(YAF_QT_H_4_MOC animworking.h mywindow.h renderwidget.h worker.h)
set(YAF_QT_UI windowbase.ui)

# povman: set more shortly path's
set(QTDIR ${YAF_USER_QT_BIN_DIR})
set(BOUNTY_ROOT ${CMAKE_SOURCE_DIR})
set(BOUNTY_INCLUDE "${CMAKE_SOURCE_DIR}/include")
set(BOUNTY_BUILD ${CMAKE_BINARY_DIR})

set(MOC_OPTIONS -DUNICODE -DWIN32 -DWIN64 -DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -D_MSC_VER=1800 -D_WIN32 -D_WIN64)

foreach(_moc_files ${YAF_QT_H_4_MOC})
    EXECUTE_PROCESS(COMMAND ${QTDIR}/moc.exe ${MOC_OPTIONS} -I${QTDIR}/mkspecs/win32-msvc2013 -I${CMAKE_SOURCE_DIR} -I${BOUNTY_ROOT} -I${BOUNTY_INCLUDE} -I${BOUNTY_BUILD} -I${QTDIR}/include -I${QTDIR}/include/QtWidgets -I${QTDIR}/include/QtGui -I${QTDIR}/include/QtCore -I. ${CMAKE_SOURCE_DIR}/src/gui/${_moc_files} -o ${CMAKE_BINARY_DIR}/moc_${_moc_files}.cpp)    
    message("MOC process file: ${CMAKE_SOURCE_DIR}/src/gui/${_moc_files}")    
    list(APPEND YAF_QT_MOC_CC ${CMAKE_BINARY_DIR}/moc_${_moc_files}.cpp)
endforeach()

EXECUTE_PROCESS(COMMAND ${QTDIR}/uic.exe ${CMAKE_SOURCE_DIR}/src/gui/windowbase.ui -o ${CMAKE_BINARY_DIR}/ui_windowbase.h)
message("UIC process file: ${CMAKE_SOURCE_DIR}/src/gui/windowbase.ui")
set(YAF_QT_UI_H ${CMAKE_BINARY_DIR}/ui_windowbase.h)

include_directories(
    "${BOUNTY_QT_DIR}/include"
    "${BOUNTY_QT_DIR}/include/QtWidgets"
    "${BOUNTY_QT_DIR}/include/QtGui"
    "${BOUNTY_QT_DIR}/include/QtCore"
    )
link_directories("${BOUNTY_QT_DIR}/lib")
#
set(QT_QTWIDGETS_LIBRARY "${BOUNTY_QT_DIR}/lib/Qt5Widgets.lib")
set(QT_QTCORE_LIBRARY    "${BOUNTY_QT_DIR}/lib/Qt5Core.lib")
set(QT_QTGUI_LIBRARY     "${BOUNTY_QT_DIR}/lib/Qt5Gui.lib")

add_library(yafarayqt SHARED ${YAF_QT_CC} ${YAF_QT_MOC_CC} ${YAF_QT_UI_H})
target_link_libraries(yafarayqt ${QT_QTWIDGETS_LIBRARY} ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} yafaraycore)
install (TARGETS yafarayqt ${YAF_TARGET_TYPE} DESTINATION ${YAF_LIB_DIR})
