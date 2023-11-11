/**
 * Block ZIP compression for BCL
 * Copyright (c) Schaudge King
 *
 * \file bzip4bcl.cpp
 */
#include <boost/cstdlib.hpp>
#include <boost/format.hpp>
#include <vector>
#include "common/Exceptions.hh"
#include "io/GzipCompressor.hh"
#include "io/FileBufWithReopen.hh"


int main(int argc, const char *argv[])
{
    uint32_t maxBlockBytes = 0xFFFF - 41;

    // 1. Read the binary bcl data into a char buffer ...
    boost::filesystem::path bclFilePath(argv[1]);
    bcl2fastq::io::FileBufWithReopen fileBuf(std::ios_base::in | std::ios_base::binary);
    fileBuf.reopen(bclFilePath, bcl2fastq::io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    int errnum = errno;
    if (!fileBuf.is_open()) {
        BOOST_THROW_EXCEPTION(bcl2fastq::common::IoError(errnum, (boost::format("Unable to open cycle BCI file '%s' for reading") % bclFilePath.string()).str()));
    }

    std::vector<char> bclBuffer;
    bclBuffer.resize(maxBlockBytes);

    std::streamsize blockSize = bcl2fastq::io::read(fileBuf, &*bclBuffer.begin(), maxBlockBytes);
    errnum = errno;
    if (blockSize != maxBlockBytes) {
        BOOST_THROW_EXCEPTION(bcl2fastq::common::InputDataError(errnum, (boost::format(
                "Unable to read header of cycle BCI file '%s': bytes_read=%d bytes_expected=%d") %
                                                                             bclFilePath.string() % blockSize %
                                                                             sizeof(maxBlockBytes)).str()));
    }
    fileBuf.close();

    /*
    // Read binary file suggested by <<Effective STL>> Item 29
    std::ifstream inputFile(argv[1], std::ios::in | std::ios::binary);
    std::vector<char> inputBuffer((std::istreambuf_iterator<char>(inputFile)), std::istreambuf_iterator<char>());
    inputFile.close();
     */

    // 2. output the bcl buffer into bgzf compressor
    std::vector<char> bclBgzfOutBuffer(bclBuffer.size());
    bcl2fastq::io::GzipCompressor compressor(bclBgzfOutBuffer,
                                             true,
                                             boost::iostreams::gzip_params(4, // config param with default=4
                                                                           boost::iostreams::zlib::deflated, // default
                                                                           15, // default
                                                                           9));

    std::streamsize bytesWritten = compressor.write(&*bclBuffer.begin(), maxBlockBytes);
    BCL2FASTQ_ASSERT_MSG( bytesWritten == maxBlockBytes,
                          "Only " << bytesWritten << " of " << bclBuffer.size() << " bytes have been written" );
    compressor.flush();

    // 3. write bcl zip block into local file
    std::ofstream outputBclBgzf(argv[2], std::ios::out | std::ios::binary);
    outputBclBgzf.write(&*bclBgzfOutBuffer.begin(), bclBgzfOutBuffer.size());
    outputBclBgzf.close();

    /*
     // write outputBuffer to .bgzf file
     std::ofstream outputFile(argv[2], std::ios::out | std::ios::bin);
     copy(outputBuffer.begin(), outputBuffer.end(), std::ostreambuf_iterator<char>(outputFile));
     outputFile.close();
     */

    return boost::exit_success;
}

