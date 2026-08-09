include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^4.1.0
)

FetchContent_Declare(
  quickstats 
  GIT_REPOSITORY https://github.com/LTLA/quickstats
  GIT_TAG master # ^0.2.0
)

FetchContent_Declare(
  sanisizer
  GIT_REPOSITORY https://github.com/LTLA/sanisizer
  GIT_TAG master # ^0.2.0
)

FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(quickstats)
FetchContent_MakeAvailable(sanisizer)
