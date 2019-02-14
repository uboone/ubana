/**
 * Tree reader input source.
 *
 * Reconstitutes art::Event products from external TTrees.
 */

#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "uboone/TreeReader/TreeReader.h"

namespace uboone {
  typedef art::Source<TreeReader> TreeReaderSource;
}

DEFINE_ART_INPUT_SOURCE(uboone::TreeReaderSource)

