###
### CMake functions to download a file
### and check its MD5 checksum
###

function(computeMD5 target md5chksum)
  execute_process(COMMAND ${CMAKE_COMMAND} -E md5sum ${target} OUTPUT_VARIABLE md5list)
  string(REGEX REPLACE "([a-z0-9]+)" "\\1;" md5list "${md5list}")
  list(GET md5list 0 md5)
  set(${md5chksum} ${md5} PARENT_SCOPE)
endfunction(computeMD5)

function(downloadAndCheckMD5 url target md5chksum)
  if (NOT ${url} EQUAL "")
    get_filename_component(targetDir ${target} PATH)
    message("Downloading ${url}...")
    file(DOWNLOAD "${url}" "${target}"
      STATUS status)
    # CMake < 2.8.10 does not seem to support HTTPS out of the box
    # and since SourceForge redirects to HTTPS, the CMake download fails
    # so we try to use Powershell (Windows) or system curl (Unix, OS X) if available
    if (NOT status EQUAL 0)
      if(WIN32)
        execute_process(COMMAND powershell -Command "(New-Object Net.WebClient).DownloadFile('${url}', '${target}')")
      else(WIN32)
        execute_process(COMMAND curl -L "${url}" -o ${target} WORKING_DIRECTORY ${targetDir})
      endif(WIN32)
    endif()
    if (NOT EXISTS ${target})
      MESSAGE(FATAL_ERROR "The download of ${url} failed.")
    endif()
    if (NOT ${md5chksum} EQUAL "")
      computeMD5(${target} md5)
      if (NOT md5 STREQUAL ${md5chksum})
        MESSAGE(FATAL_ERROR "The md5 checksum for ${target} is incorrect; expected: ${md5chksum}, found: ${md5}")
      endif()
    endif()
  endif()
endfunction(downloadAndCheckMD5)
