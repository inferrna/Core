/****************************************************************************
*      This library is free software; you can redistribute it and/or
*      modify it under the terms of the GNU Lesser General Public
*      License as published by the Free Software Foundation; either
*      version 3 of the License, or (at your option) any later version.
*
*      This library is distributed in the hope that it will be useful,
*      but WITHOUT ANY WARRANTY; without even the implied warranty of
*      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*      Lesser General Public License for more details.
*
*      You should have received a copy of the GNU Lesser General Public
*      License along with this library; if not, write to the Free Software
*      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <yafray_config.h>
#include <cstdlib>
#include <cctype>
#include <algorithm>

#ifdef WIN32
	#include <windows.h>
    #include <direct.h>
#endif

#include <core_api/scene.h>
#include <core_api/environment.h>
#include <core_api/integrator.h>
#include <core_api/imagefilm.h>
#include <yafraycore/xmlparser.h>
#include <yaf_revision.h>
#include <utilities/console_utils.h>
#include <yafraycore/imageOutput.h>

#include <gui/yafqtapi.h>

using namespace::yafaray;

int main(int argc, char *argv[])
{
	std::string xmlLoaderVersion = "TheBounty XML loader version 0.3";

	cliParser_t parse(argc, argv, 2, 1, "You need to set at least a valid XML scene file.");

	parse.setAppName(xmlLoaderVersion,
        "[OPTIONS]... <input xml file> [output filename]\n"
        "<input xml file> : A valid TheBounty XML scene file\n"
        "[output filename] : The filename of the rendered image without extension.\n"
        "*Note: If output filename is ommited the name \"bounty\" will be used instead.");

	parse.setOption("pp", "plugin-path", false,
                    "\n\tPath to load plugins."
                    "\n\t(warning: this flag overrides default path to plugins.");

	parse.setOption("vl", "verbosity-level", false,
                    "\n\tSet verbosity level, options are:\n"
                    "\t\t0 - MUTE    (Prints nothing)\n"
                    "\t\t1 - ERROR   (Prints only errors)\n"
                    "\t\t2 - WARNING (Prints only errors and warnings)\n"
                    "\t\t3 - INFO    (Prints all messages)\n");
	parse.parseCommandLine();

#ifdef RELEASE
	std::string version = std::string(VERSION);
#else
	std::string version = std::string(YAF_SVN_REV);
#endif

	renderEnvironment_t *env = new renderEnvironment_t();

	// Plugin load
	std::string ppath = parse.getOptionString("pp");
	int verbLevel = parse.getOptionInteger("vl");

	if (verbLevel >= 0) yafout.setMasterVerbosity(verbLevel);

	if (ppath.empty()) env->getPluginPath(ppath);

	if (!ppath.empty())
	{
		Y_INFO << "The plugin path is: " << ppath << yendl;
		env->loadPlugins(ppath);
	}
	else
	{
		Y_ERROR << "Getting plugin path from render environment failed!" << yendl;
		return 1;
	}

	std::vector<std::string> formats = env->listImageHandlers();

	std::string formatString = "";
	for (size_t i = 0; i < formats.size(); i++)
	{
		formatString.append("\t\t" + formats[i]);
		if (i < formats.size() - 1) formatString.append("\n");
	}

	parse.setOption("v", "version", true,
                    "\n\tDisplays this program's version.");
	parse.setOption("h", "help", true,
                    "\n\tDisplays this help text.");
	parse.setOption("op", "output-path", false,
                    "\n\tUses the path in <value> as rendered image output path.");
	parse.setOption("f", "format", false,
                    "\n\tSets the output image format, available formats are:\n\n" + formatString +
                    "\n\t\tDefault: tga.\n");
	parse.setOption("t", "threads", false,
                    "\n\tOverrides threads setting on the XML file, for auto selection use -1.");
	parse.setOption("a", "with-alpha", true,
                    "\n\tEnables saving the image with alpha channel."
                    "\n\tIf this option are omited, the 'alpha_channel' scene parameter is used.");
	parse.setOption("dp", "draw-params", true,
                    "\n\tEnables saving the image with a settings badge.");
	parse.setOption("ndp", "no-draw-params", true,
                    "\n\tDisables saving the image with a settings badge"
                    "\n\t(warning: this overrides --draw-params setting).");
	parse.setOption("cs", "custom-string", false,
                    "\n\tSets the custom string to be used on the settings badge.");
	parse.setOption("z", "z-buffer", true,
                    "\n\tEnables the rendering of the depth map (Z-Buffer)"
                    "\n\t(warning: this flag overrides XML setting).");
	parse.setOption("nz", "no-z-buffer", true,
                    "\n\tDisables the rendering of the depth map (Z-Buffer)"
                    "\n\t(warning: this flag overrides XML setting).");

	bool parseOk = parse.parseCommandLine();

	if (parse.getFlag("h"))
	{
		parse.printUsage();
		return 0;
	}

	if (parse.getFlag("v"))
	{
		Y_INFO << xmlLoaderVersion << yendl << "Built with TheBounty version " << version << yendl;
		return 0;
	}

	if (!parseOk)
	{
		parse.printError();
		parse.printUsage();
		return 0;
	}

	bool alpha = parse.getFlag("a");
	std::string format = parse.getOptionString("f");
	int threads = parse.getOptionInteger("t");
	bool drawparams = parse.getFlag("dp");
	bool nodrawparams = parse.getFlag("ndp");
	std::string customString = parse.getOptionString("cs");
	bool zbuf = parse.getFlag("z");
	bool nozbuf = parse.getFlag("nz");

	if (format.empty()) format = "tga";
	bool formatValid = false;

	for (size_t i = 0; i < formats.size(); i++)
	{
		if (formats[i].find(format) != std::string::npos) formatValid = true;
	}

	if (!formatValid)
	{
		Y_ERROR << "Couldn't find any valid image format, image handlers missing?" << yendl;
		return 1;
	}

	const std::vector<std::string> files = parse.getCleanArgs();

	if (files.size() == 0)
	{
		return 0;
	}
	// xml file name is always the first parameter without '-'
	std::string xmlFile = files[0];

    // out image file name is always the second parameter without '-'
	std::string outName = "bounty.";
    if (files.size() > 1)
    {
        outName = files[1] + "." + format;
    }
    else
    {
        size_t start_filename, filename, extension;
        std::string tmp = files[0];
        // normalize slash to UNIX style '/'
        for (int i = 0; i < tmp.length(); i++)
        {
            if (tmp[i] == '\\' || tmp[i] == '\\\\') tmp.replace(i, 1, "/");
        }
        // isolate filename from path
        if (tmp.at(tmp.length() - 1) != '/')
        {
            start_filename = tmp.rfind("/");
            tmp.erase(0, start_filename + 1);
        }        
        filename = tmp.find_last_of(".");
        extension = tmp.rfind(".");
        outName = tmp.substr(0, filename + 1) + format;
	}

	//env->Debug = debug; //disabled until proper debugging messages are set throughout the core

	// Set the full output path with filename
    // povman: this parameter need that the path are created 
    std::string outputPath = parse.getOptionString("op");

	if (outputPath.empty())
	{
		outputPath = outName;
	}
	else if (outputPath.at(outputPath.length() - 1) == '/')
	{
		outputPath += outName;
	}
	else if (outputPath.at(outputPath.length() - 1) != '/')
	{
		outputPath += "/" + outName;
	}

	scene_t *scene = new scene_t();
	env->setScene(scene);
	paraMap_t render;

	bool success = parse_xml_file(xmlFile.c_str(), scene, env, render);
	if(!success) exit(1);

	int width=320, height=240;
	int bx = 0, by = 0;
	render.getParam("width", width); // width of rendered image
	render.getParam("height", height); // height of rendered image
	render.getParam("xstart", bx); // border render x start
	render.getParam("ystart", by); // border render y start
    
    if (!alpha) render.getParam("alpha_channel", alpha);

	if(threads >= -1) render["threads"] = threads;

	if(drawparams)
	{
		render["drawParams"] = true;
		if(!customString.empty()) render["customString"] = customString;
	}

	if(nodrawparams) render["drawParams"] = false;

	if(zbuf) render["z_channel"] = true;
	if(nozbuf) render["z_channel"] = false;

	bool use_zbuf = false;
	render.getParam("z_channel", use_zbuf);

	// create output
	colorOutput_t *out = NULL;

	paraMap_t ihParams;
	ihParams["type"] = format;
	ihParams["width"] = width;
	ihParams["height"] = height;
	ihParams["alpha_channel"] = alpha;
	ihParams["z_channel"] = use_zbuf;

	imageHandler_t *ih = env->createImageHandler("outFile", ihParams);

	if(ih)
	{
		out = new imageOutput_t(ih, outputPath, bx, by);
		if(!out) return 1;
	}
	else return 1;

	if(! env->setupScene(*scene, render, *out) ) return 1;

	scene->render();
	env->clearAll();

	imageFilm_t *film = scene->getImageFilm();

	delete film;
	delete out;

	return 0;
}