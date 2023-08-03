FROM ubuntu:jammy

ARG ARCHIVE_URL="https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz"

ENV DEBIAN_FRONTEND noninteractive
ENV TEXLIVE_PATH /usr/local/texlive/current/bin/x86_64-linux

# Install dependencies
RUN apt-get update \
    && apt-get install -y \
      cpanminus build-essential wget perl \
    && apt-get clean -y \
    && apt-get autoremove -y \
    && apt-get autoclean -y

COPY ./texlive.profile ./
# Official Installation Procedure
# See: https://tug.org/texlive/quickinstall.html
RUN cd /tmp \
    && wget -nv $ARCHIVE_URL \
    && zcat < install-tl-unx.tar.gz | tar xf - \
    && cd install-tl-* \
    && perl ./install-tl --profile=texlive.profile --no-interaction

# Create a symbolic link for the installed version of texlive to a more generic path
RUN ln -s $(find /usr/local/texlive -maxdepth 1 -type d -name "20??") /usr/local/texlive/current

# Update texlive and install latex utilities
RUN export PATH=$PATH:$TEXLIVE_PATH \
    && tlmgr update --self --all \
    && tlmgr install latexindent latexmk 

# Install perl modules for latexindent
RUN cpanm Log::Log4perl Log::Dispatch::File YAML::Tiny File::HomeDir Unicode::GCString

# Add texlive path to the PATH in bash and zsh startup scripts
RUN echo "export PATH=$PATH:$TEXLIVE_PATH" >> ~/.bashrc \
    && echo "export PATH=$PATH:$TEXLIVE_PATH" >> ~/.zshrc

COPY ./.latexmkrc /root/
