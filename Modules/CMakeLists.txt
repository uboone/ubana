install(DIRECTORY . DESTINATION Modules FILES_MATCHING PATTERN "Find*.cmake")
file(COPY . DESTINATION . FILES_MATCHING PATTERN "Find*.cmake")
