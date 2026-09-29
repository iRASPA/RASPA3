#pragma once

// Which source the running binary was built from. Defined in the translation unit that
// src/raspakit/build_info.cmake writes at build time; included from the global module fragment of
// 'build_info.ixx' so that the definitions and the uses are attached to the same (global) module.
// Read these through the 'build_info' module (BuildInfo::version() etc.), not directly.
namespace BuildInfoData
{
extern const char *const version;      // project version, e.g. "3.1.0"
extern const char *const commit;       // abbreviated git hash, or "unknown"
extern const char *const commitDate;   // ISO date of that commit, or "unknown"
extern const char *const treeState;    // "" or ", with uncommitted changes"
extern const char *const compileDate;  // __DATE__ __TIME__ of the last relink-worthy change
}  // namespace BuildInfoData
