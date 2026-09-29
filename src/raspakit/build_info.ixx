module;

#include "build_info.h"

export module build_info;

import std;

/**
 * \namespace BuildInfo
 * \brief Which source the running binary was built from.
 *
 * The values come from the translation unit that 'build_info.cmake' regenerates at build time from the git
 * state of the source tree, so a result can be traced back to the commit that produced it. The same text is
 * shown by 'raspa3 --version' and written at the top of every output file.
 */
export namespace BuildInfo
{
/// Project version, e.g. "3.1.0".
inline std::string_view version() { return BuildInfoData::version; }

/// Abbreviated git hash of the commit the binary was built from, or "unknown" outside a git checkout.
inline std::string_view commit() { return BuildInfoData::commit; }

/// ISO date of that commit, or "unknown".
inline std::string_view commitDate() { return BuildInfoData::commitDate; }

/// Empty, or ", with uncommitted changes" when tracked files differed from the commit at build time.
inline std::string_view treeState() { return BuildInfoData::treeState; }

/// __DATE__ __TIME__ of the last build in which the commit or the dirty state changed.
inline std::string_view compileDate() { return BuildInfoData::compileDate; }

/// The multi-line summary shown by 'raspa3 --version' and 'raspa3 --help'.
inline std::string summary()
{
  return std::format("raspa3 {}\n  commit:   {} ({}){}\n  compiled: {}", version(), commit(), commitDate(),
                     treeState(), compileDate());
}
}  // namespace BuildInfo
