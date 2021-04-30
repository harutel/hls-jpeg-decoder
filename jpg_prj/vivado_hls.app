<project xmlns="com.autoesl.autopilot.project" name="jpg_prj" top="JpegDecodeHW">
    <files>
        <file name="src/loadjpg.cpp" sc="0" tb="false" cflags=""/>
        <file name="../../src/main.cpp" sc="0" tb="1" cflags=" "/>
        <file name="../../src/openjpg.cpp" sc="0" tb="1" cflags=""/>
    </files>
    <includePaths/>
    <libraryPaths/>
    <Simulation>
        <SimFlow name="csim" clean="true" csimMode="0" lastCsimMode="0" compiler="true"/>
    </Simulation>
    <solutions xmlns="">
        <solution name="01_no_opt" status="inactive"/>
        <solution name="99_winnder" status="active"/>
    </solutions>
</project>

