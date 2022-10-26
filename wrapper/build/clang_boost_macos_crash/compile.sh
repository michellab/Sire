

clang++ -DBOOST_PYTHON_DYN_LINK -DBOOST_PYTHON_NO_LIB -I$CONDA_PREFIX/Library/include -I$CONDA_PREFIX/include/python3.7m -isystem $CONDA_PREFIX/include $CONDA_PREFIX/./mkspecs/macx-clang -Os -DNDEBUG -Wall -Wno-attributes -pipe -w -undefined dynamic_lookup -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk -mmacosx-version-min=10.9 -fPIC -fPIC -MD -MT test_bp.cpp.o -MF test_bp.cpp.o.d -o test_bp.cpp.o -c test_bp.cpp

clang++ -DNDEBUG -Wall -Wno-attributes -pipe -w -undefined dynamic_lookup -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk -mmacosx-version-min=10.9 -dynamiclib -Wl,-headerpad_max_install_names -compatibility_version 2022.0.0 -current_version 2022.2.0 -o _test_bp.so test_bp.cpp.o -L$CONDA_PREFIX/lib -lboost_python37
