#!/usr/bin/env python3
import sys, os, tempfile, shutil, gzip, selectors, re

# If called directly, this works fine
try:
    from collapse import collapse
    from clip import clip
    from trim import trim
    from configutator import ConfigMap, ArgMap, loadConfig
    from valitator import Path, PathOrNone, Executable
# If installed
except ImportError:
    from ProDuSe.collapse import collapse
    from ProDuSe.clip import clip
    from ProDuSe.trim import trim
    from ProDuSe.configutator import ConfigMap, ArgMap, loadConfig
    from ProDuSe.valitator import Path, PathOrNone, Executable

from multiprocessing import Process
from subprocess import Popen, PIPE


def pipeOpener(file, flags):
    return os.open(file, flags|os.O_RDONLY|os.O_NONBLOCK)

stripDecorators = re.compile("\x1b[^A-HJKSTfminsulh]*.|\t|\r|\n")

def displayVerbose(screen, sel: selectors.DefaultSelector, data, poll):
    from asciimatics.screen import Screen, StopApplication, NextScene, ResizeScreenError
    from asciimatics.scene import Scene
    import asciimatics.widgets as Widgets
    from asciimatics.effects import Matrix, _Trail, randint
    from time import time

    # It was too good to not include
    class _BioTrail(_Trail):
        def update(self, reseed):
            chars = ['A', 'G', 'C', 'T']
            if self._clear:
                for i in range(0, 3):
                    self._screen.print_at(" ",
                                          self._x,
                                          self._screen.start_line + self._y + i)
                self._maybe_reseed(reseed)
            else:
                for i in range(0, 3):
                    self._screen.print_at(chars[randint(0,3)],
                                          self._x,
                                          self._screen.start_line + self._y + i,
                                          Screen.COLOUR_GREEN)
                for i in range(4, 6):
                    self._screen.print_at(chars[randint(0,3)],
                                          self._x,
                                          self._screen.start_line + self._y + i,
                                          Screen.COLOUR_GREEN,
                                          Screen.A_BOLD)
                self._maybe_reseed(reseed)

    class BioMatrix(Matrix):
        def reset(self):
            self._chars = [_BioTrail(self._screen, x) for x in
                           range(self._screen.width)]

    Widgets.Frame.palette['section_header'] = (Screen.COLOUR_GREEN, Screen.A_UNDERLINE, Screen.COLOUR_BLUE)

    window = Widgets.Frame(screen, 30, 200, title="ProDuSe V2.0", data=data)
    scene = Scene([BioMatrix(screen), window])
    layout = Widgets.Layout([1])
    window.add_layout(layout)
    outputFile = Widgets.Text('Output', 'output')
    outputFile.disabled = True
    layout.add_widget(outputFile)
    # Trim
    label = Widgets.Label("Trim:")
    label.custom_colour = 'section_header'
    layout.add_widget(label)
    trim = Widgets.TextBox(1, name='trim')
    trim.disabled = True
    layout.add_widget(trim)
    layout.add_widget(Widgets.Divider(False))

    #BWA
    label = Widgets.Label("BWA:")
    label.custom_colour = 'section_header'
    layout.add_widget(label)
    box = Widgets.TextBox(len(data['bwa']), name='bwa')
    box.disabled = True
    layout.add_widget(box)
    layout.add_widget(Widgets.Divider(False))

    # Clip
    label = Widgets.Label("Clip:")
    label.custom_colour = 'section_header'
    layout.add_widget(label)
    box = Widgets.TextBox(1, name='clip')
    box.disabled = True
    layout.add_widget(box)
    layout.add_widget(Widgets.Divider(False))

    # Collapse
    label = Widgets.Label("Collapse:")
    label.custom_colour = 'section_header'
    layout.add_widget(label)
    box = Widgets.TextBox(1, name='collapse')
    box.disabled = True
    layout.add_widget(box)
    layout.add_widget(Widgets.Divider(False))

    window.fix()
    screen.set_scenes([scene])
    screen.draw_next_frame()
    lastTime = time()
    while poll():
        for key, mask in sel.select(0.5):
            key.data.pop(0)
            key.data.append(stripDecorators.sub(lambda m: '   ' if m.group(0) == '\t' else '', key.fileobj.readline().decode()))
        if time() - lastTime > 1/20: # Dump output until a reasonable amount of time has passed to avoid the buffer filling and slowing down everything
            if screen.has_resized():
                raise ResizeScreenError('')
            window.data = data
            screen.force_update()
            screen.draw_next_frame()
            lastTime = time()

    raise StopApplication('')

@ConfigMap(fastq1='fastqs[0]', fastq2='fastqs[1]', output='join(``, [output, name, `.bam`])', verbose=None)
@ArgMap(_func=None)
def ProDuSe(fastq1: Path, fastq2: Path, reference: Path, output: str, bwa: Executable, samtools: Executable, bed: PathOrNone = None, verbose: bool=False):
    """
    Manages the subprocesses and interprocess communication of Trim, BWA, Clip, samtools, Collapse, and Filter
    :param fastq1: Path to the first fastq
    :param fastq2: Path to the second fastq
    :param reference: Path to the reference fasta
    :param bed: Path to a bed file containing regions to restrict reads to (optional)
    :param output: Path to output file
    :param bwa: Path to bwa executable and any added parameters ex. mem or -t
    :param samtools: Path to samtools executable
    :param verbose: Provide verbose output while processing
    :return: None
    """
    global config
    stderr = sys.stderr
    temp_dir = tempfile.mkdtemp()
    if verbose: stderr.write("Pipes created in {}\n".format(temp_dir))
    trimOutPath = os.path.join(temp_dir, 'trimOut')
    os.mkfifo(trimOutPath)

    debugPipes = []

    # Check if fastqs gzipped
    inFile1 = openGZ(fastq1, 'rb')
    inFile2 = openGZ(fastq2, 'rb')

    # Start bwa subprocess
    if verbose: stderr.write("Starting BWA subprocess..\n")
    bwa = Popen(bwa.split(' ') + ['-Cp', reference, trimOutPath], stdout=PIPE, stderr=PIPE if verbose else None)
    clipIn = bwa.stdout

    # Start trim subprocesses
    if verbose:
        stderr.write("Starting Trim subprocesses..\n")
        r, w = os.pipe()
        os.set_inheritable(w, True)
        debugPipes.append(open(r, 'rb'))
    trimOut = open(trimOutPath, 'wb')
    config[trim].update(dict(inStream=inFile1, mateStream=inFile2, outStream=trimOut, logStream=open(w if verbose else os.devnull, 'w', buffering=1)))
    trimProc = Process(target=trim, kwargs=config[trim])
    trimProc.start()
    trimOut.close()

    debugPipes.append(bwa.stderr)

    if bed:
        #Start samtools view to restrict to coordinates
        if verbose: stderr.write("Starting samtools view subprocess..\n")
        view = Popen([samtools, "view", "-uL", bed, "-"], stdin=clipIn, stdout=PIPE)
        clipIn.close()
        clipIn = view.stdout

    # Start sort by coord
    if verbose: stderr.write("Starting sort subprocess..\n")
    sort = Popen([samtools, "sort", "-l0", "-"], stdout=PIPE, stdin=PIPE)

    #Start clipping subprocess
    if verbose:
        stderr.write("Starting clipping subprocess..\n")
        r, w = os.pipe()
        os.set_inheritable(w, True)
        debugPipes.append(open(r, 'rb'))
    clipProc = Process(target=clip, args=(clipIn, sort.stdin), kwargs=dict(alternate=True, verbose=verbose, logStream=open(w if verbose else os.devnull, 'w', buffering=1)))
    clipProc.start()
    clipIn.close()
    sort.stdin.close()

    if verbose:
        stderr.write("Starting collapse subprocess..\n")
        r, w = os.pipe()
        os.set_inheritable(w, True)
        debugPipes.append(open(r, 'rb'))
    config[collapse].update(dict(inStream=sort.stdout, outStream=open(output, 'wb+'), logStream=open(w if verbose else os.devnull, 'w', buffering=1)))
    collapseProc = Process(target=collapse, kwargs=config[collapse])
    collapseProc.start()
    sort.stdout.close()

    if verbose:
        from collections import OrderedDict
        sel = selectors.DefaultSelector()
        lineBuffers = [('trim', ['']), ('bwa', ['' for _ in range(10)]), ('clip', ['']), ('collapse', [''])]
        for i in range(len(debugPipes)):
            sel.register(debugPipes[i], selectors.EVENT_READ, lineBuffers[i][1])
        poll = lambda : trimProc.is_alive() or clipProc.is_alive() or bwa.poll() or sort.poll() or collapseProc.is_alive()
        from asciimatics.screen import Screen
        from asciimatics.exceptions import ResizeScreenError, StopApplication
        # Handle screen resize by rebuilding screen
        while True:
            try:
                Screen.wrapper(displayVerbose, arguments=(sel, OrderedDict(lineBuffers + [('output', output)]), poll))
            except ResizeScreenError:
                continue
            except StopApplication:
                break

    trimProc.join()
    inFile1.close()
    inFile2.close()
    clipProc.join()
    collapseProc.join()
    bwa.wait()
    if bed:
        view.wait()
    sort.wait()
    shutil.rmtree(temp_dir)

<<<<<<< HEAD:ProDuSe/produse.py
def main(args=None):
    ConfigMap(_func='trim', verbose=None, logStream=None)(trim)
    ConfigMap(_func='collapse', verbose=None, logStream=None)(collapse)
    ArgMap(_func=trim.__name__, verbose='verbose', logStream=None)(trim)
    ArgMap(_func=collapse.__name__, verbose='verbose', logStream=None)(collapse)
    if args is None:
        args = sys.argv
    cfgs = loadConfig(sys.argv, (ProDuSe, trim, collapse), batchExpression='samples')
    global config
    while True:
        try:
            argmap = next(cfgs)
            config = argmap
            ProDuSe(**argmap[ProDuSe])
        except ValueError:
            continue
        except StopIteration:
            break

if __name__ == "__main__":
    main()

=======

def main(args=None):

    if args is None:
        args = sys.argv
    ConfigMap(_func='trim', verbose=None)(trim)
    ConfigMap(_func='collapse', verbose=None)(collapse)
    ArgMap(_func=trim.__name__, verbose='verbose')(trim)
    ArgMap(_func=collapse.__name__, verbose='verbose')(collapse)
    for argmap, params in loadConfig(args, (ProDuSe, trim, collapse), batchExpression='samples'):
        global config
        config = argmap
        ProDuSe(**argmap[ProDuSe])

if __name__ == "__main__":
    main()
>>>>>>> c1e0b8f7405451bb9e981e9490742ed78e82c2b3:ProDuSe/produse.py
