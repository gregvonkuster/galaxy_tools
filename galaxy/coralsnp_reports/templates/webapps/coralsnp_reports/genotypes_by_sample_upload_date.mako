<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">Genotypes of all uploaded samples</h3>
        <h4 align="center">
            Click the number of genotypes to view the number of genotypes for uploaded samples per month
        </h4>
        <table align="center" class="colored">
            %if num_genotypes == 0:
                <tr><td>There are no genotypes</td></tr>
            %else:
                <tr class="header">
                    <td align="center">
                        Number of genotypes for all uploaded samples
                    </td>
                </tr>
                <tr class="tr">
                    <td align="center">
                        <a href="${h.url_for( controller='genotypes', action='per_month', sort_id='default', order='default' )}">
                            ${num_genotypes}
                        </a>
                    </td>
                </tr>
            %endif
        </table>
    </div>
</div>

