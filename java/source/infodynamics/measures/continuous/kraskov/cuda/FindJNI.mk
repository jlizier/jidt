ifeq ($(JAVA_HOME),)
JAVAC ?= $(realpath $(call which,javac))
JAVA_HOME = $(abspath $(dir $(JAVAC))..)
endif

ifneq ($(JAVA_HOME),)
JNI_INCLUDE ?= $(JAVA_HOME)/include
endif

ifeq ($(JNI_INCLUDE),)
$(error could not determine JNI include dir, try specifying either \
    JAVA_HOME or JNI_INCLUDE)
endif

TARGETTRIPLET := $(shell $(CC) -dumpmachine)
ifeq ($(findstring mingw,$(TARGETTRIPLET)),mingw)
JNI_PLATFORM:= win32
else
ifeq ($(findstring linux,$(TARGETTRIPLET)),linux)
JNI_PLATFORM:= linux
endif
endif

JNI_PLATFORM_INCLUDE ?= $(JNI_INCLUDE)/$(JNI_PLATFORM)

