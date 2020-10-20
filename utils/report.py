import os


def str2fileBaseName(strName):
    return '_'.join(strName.split(' '))

def get_rows(flatList, nCol):
    rowList = []
    if len(flatList) < nCol:
        rowList = [flatList]
    else:
        currRow = []
        for i in range(len(flatList)):
            currRow.append(flatList[i])
            if len(currRow) == nCol or\
                    i == len(flatList)-1:
                rowList.append(currRow)
                currRow = []
    return rowList

def tex_header(projName, author):
    tmp = '\\documentclass[12pt]{article}\n'
    tmp+= '\\usepackage{graphicx}\n'
    tmp+= '\\usepackage[margin=0.2in]{geometry}\n'
    tmp+= '\\usepackage[section]{placeins}\n'
    tmp+= '\\usepackage{morefloats}\n'
    tmp+= '\\title{%s}\n'%projName
    tmp+= '\\author{%s}\n'%author
    tmp+= '\\date{\\today}\n'
    tmp+= '\\begin{document}\n\\maketitle\n\\FloatBarrier\n'
    return tmp

def tex_end():
    return '\\end{document}'

def tex_IMG(imgName, width):
    print('write '+imgName)
    tmp = '\\begin{minipage}{%f\\textwidth}\n\\centering\n'%width
    tmp+= '\\includegraphics[width=.99\\linewidth]{%s}\n\\end{minipage}\n'%imgName
    return tmp

def tex_IMG_row(imgRow, nCol):
    tmp = '\\begin{figure}[!htb]\n'
    for img in imgRow:
        tmp+= tex_IMG(img, 0.95/nCol)
    tmp += '\\end{figure}\n'
    tmp += '\\FloatBarrier\n'
    return tmp

def tex_space():
    return '\\FloatBarrier\n'

def tex_newPage():
    return '\\newpage\n'

def tex_splitLine():
    return '\\rule{\\textwidth}{0.2mm}\n'

def render_pdf(baseName, tex, rmTemp=True):
    f = open(baseName+'.tex', 'w')
    f.write(tex)
    f.close()
    os.system('pdflatex %s.tex'%baseName)
    if rmTemp:
        os.system('rm %s.log %s.aux'%(baseName, baseName))

def gen_pdf_test(
    projName = 'Calculation Result',
    author = 'Zisheng Zhang',
    description = '',
    imgList = None,
    nCol = 4,
    rmTemp=True,
    ):
    story = tex_header(projName, author)
    story+= '%s\n'%description
    imgRow = get_rows(imgList, nCol)
    for r in imgRow:
        story+= tex_IMG_row(r, nCol=nCol)
    story+= tex_end()
    baseName = str2fileBaseName(projName)
    render_pdf(baseName, story, rmTemp)

def gen_pdf_test2(
    projName = 'Calculation Result',
    author = 'Zisheng Zhang',
    description = '',
    bigFig = None,
    imgList = None,
    nCol = 4,
    rmTemp=True,
    ):
    story = tex_header(projName, author)
    story+= '%s\n'%description
    story+= tex_splitLine()
    for i in bigFig:
        story+= tex_IMG(i, 1.0)
        story += tex_space()
    story+= tex_newPage()
    imgRow = get_rows(imgList, nCol)
    for r in imgRow:
        story+= tex_IMG_row(r, nCol=nCol)
    story+= tex_end()
    baseName = str2fileBaseName(projName)
    render_pdf(baseName, story, rmTemp)